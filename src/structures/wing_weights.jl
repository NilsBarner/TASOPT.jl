"""
    size_wing_section!(section, sweep, sigfac)

Calculates Loads and thicknesses for wing sections

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `section::TASOPT.structures.Wing.WingSection`: Wing Section to be sized.
    - `sweep::Float64`: Wing sweep.
    - `sigfac::Float64`: Stress factor
"""
function size_wing_section!(section, sweep, sigfac)
    shear_load = section.max_shear_load
    moment = section.moment

    cross_section = section.cross_section

    Eweb = section.webs.material.E
    Ecap = section.caps.material.E
    Gcap = section.caps.material.G
    Gweb = section.webs.material.G
    sigcap = section.caps.material.σmax * sigfac
    tauweb = section.webs.material.τmax * sigfac

    cs = section.co
    cosL = cosd(sweep)
    cp = cs*cosL

    h_avg, h_rms = get_average_sparbox_heights(section.cross_section)

    web_height = cross_section.web_to_box_height * cross_section.thickness_to_chord

    tbweb, Abweb = size_web(tauweb, shear_load, cp, web_height)
    
    tbcap, Abcap = size_cap(sigcap, moment, cross_section.thickness_to_chord,
        cross_section.width_to_chord, h_rms, cs, cosL)

    # EI_xx
    section.EI[1] = Ecap * cp^4 * (h_rms^3 - (h_rms - 2.0 * tbcap)^3) * cross_section.width_to_chord / 12.0 +
                                 Eweb * cp^4 * tbweb * web_height^3 / 6.0
    # EI_yy 
    section.EI[4] = Ecap * cp^4 * tbcap * cross_section.width_to_chord^3 / 6.0 +
                                Eweb * cp^4 * tbweb * web_height * 0.5 * cross_section.width_to_chord^2
    # println("GJ SURFW: cp = $cp, width_to_chord = $(cross_section.width_to_chord), tbweb = $tbweb, h_avg = $h_avg, tbcap = $tbcap, web_to_box_height = $(cross_section.web_to_box_height), thickness_to_chord = $(cross_section.thickness_to_chord), Gweb = $Gweb, Gcap = $Gcap")
    section.GJ = cp^4 * 2.0*((cross_section.width_to_chord-tbweb)*(h_avg-tbcap))^2 /
        (  (cross_section.web_to_box_height*section.cross_section.thickness_to_chord-tbcap)/(Gweb*tbweb) +
        (   cross_section.width_to_chord -tbweb)/(Gcap*tbcap) )

    return tbweb, tbcap, Abcap, Abweb  
end

"""
    wing_weights!(wing, po, gammat, gammas, 
       Nload, We, neout, dyeout, neinn, dyeinn, sigfac, rhofuel; n_wings=2.0)

Calculates Wing or Tail loads, stresses, weights of individual wing sections.
Also returns the material gauges, torsional and bending stiffness. Formerly, `get_wing_weights!()`.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `wing::TASOPT.structures.Wing`: Wing structure.
    - `po::Float64`: Wing root loading magnitude.
    - `gammat::Float64`: Wing outer section taper ratio.
    - `gammas::Float64`: Wing inner section taper ratio.
    - `Nload::Int`: Max vertical load factor for wing bending loads.
    - `We::Float64`: Weight of the engine.
    - `neout::Int`:  Number of engines mounted outboard of the wing break (strut attachment point).
    - `dyeout::Float64`: Spanwise moment arm for outboard engines, measured from the wing break station.
    - `neinn::Int`: Number of engines mounted inboard of the wing break (strut attachment point).
    - `dyeinn::Float64`: Spanwise moment arm for inboard engines, measured from the wing break station.
    - `sigfac::Float64`: Stress Factor.
    - `rhofuel::Float64`: Density of the fuel.
    - `n_wings::Int64`: Number of total wings (1 for Vtail).

See [Geometry](@ref geometry),  [Wing/Tail Structures](@ref wingtail), and Section 2.7  of the [TASOPT Technical Description](@ref dreladocs). 
"""
function wing_weights!(wing, po, gammat, gammas,
       Nload, We, neout, dyeout, neinn, dyeinn, sigfac, rhofuel; n_wings=2.0)

    Wpointload = sum(-point_load.force[3] for point_load in wing.point_loads)  # NILS: modified to sum over all point loads ("-" sign accounts for positive z-direction being upwards, whereas weight is positive when force pointing downwards)

    # NILS' logic for adding wing point load
    ywingbreak = wing.layout.ηs * wing.layout.span / 2.0
    Wpointloadout = 0.0
    Wpointloadin = 0.0
    dyWpointloadout = 0.0
    dyWpointloadin = 0.0
    for point_load in wing.point_loads
        Wpointload = -point_load.force[3]
        if point_load.r[2] isa Number
            ypointload = point_load.r[2] * wing.layout.span / 2.0
        elseif point_load.r[2] isa Dict
            ypointload = point_load.r[2]["frac_span"] * wing.layout.span / 2.0
        end
        if ypointload > ywingbreak
            Wpointloadout += Wpointload
            dyWpointloadout += Wpointload * (ypointload - ywingbreak)
        elseif ypointload <= ywingbreak
            Wpointloadin += Wpointload
            dyWpointloadin += Wpointload * (ypointload - ywingbreak)  # SIGN TO BE CHECKED
        end
    end

    # println("Wpointloadin, Wpointloadout = ", Wpointloadin, Wpointloadout)

    # NILS: place-holders (assume zero x-moment from point loads for now - TO BE INCLUDED)
    dxWpointloadout = 0.0
    dxWpointloadin = 0.0

    tauweb,sigstrut = wing.inboard.webs.material.τmax * sigfac, wing.strut.material.σmax * sigfac

    cosL = cosd(wing.layout.sweep)
    sinL = sind(wing.layout.sweep)

    # Calculate non-dim span coordinate at span break and root (ηs and ηo resp.)
    etao = wing.layout.ηo
    etas = wing.layout.ηs

    # Tip roll off Lift (modeled as a point load) and it's moment about ηs
    dLt = wing.tip_lift_loss * po * wing.layout.root_chord * gammat * wing.outboard.λ
    dMt = dLt * 0.5 * wing.layout.span * (1.0 - etas)

    h_avgo, h_rmso = get_average_sparbox_heights(wing.inboard.cross_section)
    h_avgs, h_rmss = get_average_sparbox_heights(wing.outboard.cross_section)

    # Outboard section:
    #---- strut-attach shear,moment from outer-wing loading. 
    #     Note added term to account for any outboard engines.
    #     If neout = 0 this simplifies to Drela's version which assumes engine
    #     fixed at ηs locations only.
    wing.outboard.max_shear_load = (po * wing.layout.span / 4.0) * (gammas + gammat) * (1.0 - etas) +
                                   dLt - Nload * wing.outboard.weight - Nload * neout * We -
                                   Nload * Wpointloadout  # NILS
    wing.outboard.moment = (po * wing.layout.span^2 / 24.0) * (gammas + 2.0 * gammat) * (1.0 - etas)^2 +
                           dMt - Nload * wing.outboard.dyW - Nload * neout * We * dyeout -
                           Nload * dyWpointloadout  # NILS

    # println("Nload * neout * We, Nload * Wpointloadout = ", Nload * neout * We, " ", Nload * Wpointloadout)
    # println("Nload * neout * We * dyeout, Nload * dyWpointloadout = ", Nload * neout * We * dyeout, " ", Nload * dyWpointloadout)

    #---- size strut-attach station at etas
    cs = wing.layout.root_chord*wing.inboard.λ

    tbwebs, tbcaps, Abcaps, Abwebs = size_wing_section!(wing.outboard, wing.sweep, sigfac)
    # Inboard Section:
    if(!wing.has_strut) 
        #----- no strut, with or without engine at etas
        ls = 0.
        Tstrut = 0.
        Rstrut = 0.
        Pstrut = 0.

        # Modifed to account for bending relief from multiple engines.
        # dyeinn allows engine to be in locations other than ηs
        So = wing.outboard.max_shear_load - Nload*neinn*We +
            0.25*po*wing.layout.span*(1.0+gammas)*(etas-etao) -
            Nload*wing.inboard.weight -
            Nload * Wpointloadin  # NILS
        # println("Nload*neinn*We, Nload * Wpointloadin = ", Nload*neinn*We, " ", Nload * Wpointloadin)
        Mo = wing.outboard.moment + wing.outboard.max_shear_load*0.5*wing.layout.span*(etas-etao) +
            (1.0/24.0)*po*wing.layout.span^2*(1.0+2.0*gammas)*(etas-etao)^2 -
            Nload*wing.inboard.dyW - Nload*neinn*We*dyeinn -
            Nload * dyWpointloadin  # NILS
        # println("Nload*neinn*We*dyeinn, Nload * dyWpointloadin = ", Nload*neinn*We*dyeinn, " ", Nload * dyWpointloadin)
        # println()

        #----- limit So,Mo to Ss,Ms, which might be needed with heavy outboard engine
        #-      (rules out negatively-tapered structure, deemed not feasible for downloads)
        wing.inboard.max_shear_load = max(So ,wing.outboard.max_shear_load)
        wing.inboard.moment = max(Mo ,wing.outboard.moment)

        tbwebo, tbcapo, Abcapo, Abwebo = size_wing_section!(wing.inboard, wing.sweep, sigfac)

        lsp = 0.
        Tstrutp = 0.
        wing.strut.cos_lambda = 1.0

    else
        #----- strut present
        ls = sqrt(wing.strut.z^2 + (0.5*wing.layout.span*(etas-etao))^2)
        Rstrut = (po*wing.layout.span/12.0)*(etas-etao)*(1.0+2.0*gammas) + wing.outboard.max_shear_load
        Tstrut = Rstrut*ls/wing.strut.z
        #c     Pstrut = Rstrut*0.5*wing.layout.span*(etas-etao)/zs

        #----- inboard shear,moment used for sparbox sizing
        wing.inboard.max_shear_load = wing.outboard.max_shear_load
        wing.inboard.moment = wing.outboard.moment
        #
        #----- size inboard station at etao
        tbwebo, tbcapo, Abcapo, Abwebo = size_wing_section!(wing.inboard, wing.sweep, sigfac)

        #----- total strut length, tension
        lsp = sqrt(wing.strut.z^2 + (0.5*wing.layout.span*(etas-etao)/cosL)^2)
        Tstrutp = Tstrut*lsp/ls
        wing.strut.cos_lambda = ls/lsp
        wing.strut.axial_force = Tstrutp/sigstrut
        wing.strut.weight   = 2.0*wing.strut.material.ρ*gee*wing.strut.axial_force*lsp
        wing.strut.dxW = wing.strut.weight * 0.25*wing.layout.span*(etas-etao) * sinL/cosL
    end

    Abfuels = calc_sparbox_internal_area(wing.inboard.cross_section.width_to_chord, h_avgs, tbcaps, tbwebs)
    Abfuelo = calc_sparbox_internal_area(wing.inboard.cross_section.width_to_chord, h_avgo, tbcapo, tbwebo) 
    
    

    Vcen = wing.layout.root_chord^2*wing.layout.span*  etao / 2.0

    Vinn = wing.layout.root_chord^2*wing.layout.span* (etas-etao) *
        (1.0 + wing.inboard.λ + wing.inboard.λ^2)/6.0 *
        cosL
    Vout = wing.layout.root_chord^2*wing.layout.span* (1.0 -etas) *
        (wing.inboard.λ^2 + wing.inboard.λ*wing.outboard.λ + wing.outboard.λ^2)/6.0 *
        cosL

    dxVinn = wing.layout.root_chord^2*wing.layout.span^2 * (etas-etao)^2 *
        (1.0 + 2.0*wing.inboard.λ + 3.0*wing.inboard.λ^2)/48.0 *
        sinL
    dxVout = wing.layout.root_chord^2*wing.layout.span^2 * (1.0 -etas)^2 *
        (wing.inboard.λ^2 + 2.0*wing.inboard.λ*wing.outboard.λ + 3.0*wing.outboard.λ^2)/48.0 *
        sinL +
        wing.layout.root_chord^2*wing.layout.span^2 * (etas-etao)*(1.0 -etas) *
        (wing.inboard.λ^2 + wing.inboard.λ*wing.outboard.λ + wing.outboard.λ^2)/12.0 *
        sinL

    dyVinn = wing.layout.root_chord^2*wing.layout.span^2 * (etas-etao)^2 *
        (1.0 + 2.0*wing.inboard.λ + 3.0*wing.inboard.λ^2)/48.0 *
        cosL
    dyVout = wing.layout.root_chord^2*wing.layout.span^2 * (1.0 -etas)^2 *
        (wing.inboard.λ^2 + 2.0*wing.inboard.λ*wing.outboard.λ + 3.0*wing.outboard.λ^2)/48.0 *
        cosL

    #---- set chord^2 weighted average areas for inner panel
    Abcapi = (Abcapo + Abcaps*wing.inboard.λ^2)/(1.0+wing.inboard.λ^2)
    Abwebi = (Abwebo + Abwebs*wing.inboard.λ^2)/(1.0+wing.inboard.λ^2)

    Wscen   = (wing.inboard.caps.ρ*Abcapo + wing.inboard.webs.ρ*Abwebo)*gee*Vcen
    Wsinn   = (wing.inboard.caps.ρ*Abcapi + wing.inboard.webs.ρ*Abwebi)*gee*Vinn
    Wsout   = (wing.inboard.caps.ρ*Abcaps + wing.inboard.webs.ρ*Abwebs)*gee*Vout

    dxWsinn = (wing.inboard.caps.ρ*Abcapi + wing.inboard.webs.ρ*Abwebi)*gee*dxVinn
    dxWsout = (wing.inboard.caps.ρ*Abcaps + wing.inboard.webs.ρ*Abwebs)*gee*dxVout

    dyWsinn = (wing.inboard.caps.ρ*Abcapi + wing.inboard.webs.ρ*Abwebi)*gee*dyVinn + dyWpointloadin  # NILS: added last term (placed here as gets returned separately per section - compare with dxWsinn)
    dyWsout = (wing.inboard.caps.ρ*Abcaps + wing.inboard.webs.ρ*Abwebs)*gee*dyVout + dyWpointloadout  # NILS: added first term (placed here as gets returned separately per section - compare with dxWsout)


    Abfueli= (Abfuelo + Abfuels*wing.inboard.λ^2)/(1.0+wing.inboard.λ^2)

    Wfcen   = rhofuel*Abfuelo *gee*Vcen
    Wfinn   = rhofuel*Abfueli *gee*Vinn
    Wfout   = rhofuel*Abfuels *gee*Vout

    # NILS: calculate wing box volumes like
    # in fuel weight calculations above
    # NOTE: while the terms match those in
    # equation (570) of tasopt.pdf, it seems
    # odd that Abfuelo gets multiplied by Vcen
    # and Abfuels by Vout -> perhaps a BUG???
    if isdefined(wing, :center)
        wing.center.volume = Abfuelo * Vcen  # NILS: undefined during first iteration
    end
    wing.inboard.volume = Abfueli * Vinn
    wing.outboard.volume = Abfuels * Vout

    dxWfinn = rhofuel*Abfueli *gee*dxVinn
    dxWfout = rhofuel*Abfuels *gee*dxVout

    dyWfinn = rhofuel*Abfueli *gee*dyVinn
    dyWfout = rhofuel*Abfuels *gee*dyVout

    wing.inboard.caps.weight.W   = 2.0*wing.inboard.caps.ρ*gee*(Abcapo*Vcen + Abcapi*Vinn + Abcaps*Vout)
    wing.inboard.webs.weight.W   = 2.0*wing.inboard.webs.ρ*gee*(Abwebo*Vcen + Abwebi*Vinn + Abwebs*Vout)

    dxWcap = 2.0*wing.inboard.caps.ρ*gee*( Abcapi*dxVinn + Abcaps*dxVout )
    dxWweb = 2.0*wing.inboard.webs.ρ*gee*( Abwebi*dxVinn + Abwebs*dxVout )

    Vout = wing.layout.root_chord^2*wing.layout.span* (1.0 -etas) *
    (wing.inboard.λ^2 + wing.inboard.λ*wing.outboard.λ + wing.outboard.λ^2)/6.0 *
    cosL

    wing.inboard.caps.thickness = tbcapo
    wing.inboard.webs.thickness = tbwebo
    wing.outboard.caps.thickness = tbcaps
    wing.outboard.webs.thickness = tbwebs

    fwadd = wing_additional_weight(wing)
    Wwing = n_wings * (Wscen + Wsinn + Wsout) * (1.0 + fwadd) + Wpointloadin + Wpointloadout  # NILS: added last two terms (placed here rather than with individual sections to avoid applying `* (1.0 + fwadd)` to point loads)
    wing.dxW = n_wings * (dxWsinn + dxWsout) * (1.0 + fwadd) + dxWpointloadin + dxWpointloadout  # NILS: added last two terms (placed here rather than with individual sections to avoid applying `* (1.0 + fwadd)` to point loads)
    # println("Wpointloadin + Wpointloadout = ", Wpointloadin + Wpointloadout)

    return Wwing, Wscen, Wsinn, Wsout, dyWsinn, dyWsout, Wfcen, Wfinn, Wfout,  # Wscen added by NILS for plotting
    dxWfinn, dxWfout, dyWfinn, dyWfout, lsp

end # wing_weights


"""
    size_cap(σmax, moment, h̄, w̄, h_rms, c, cosΛ)

Calculates area and thickness of wing caps based on maximum stress and geometric parameters.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `σmax::Float64`: Maximum allowable stress
    - `moment::Float64`: Bending moment
    - `h̄::Float64`: Height-to-chord ratio
    - `w̄::Float64`: Width-to-chord ratio
    - `h_rms::Float64`: RMS height of sparbox
    - `c::Float64`: Chord length
    - `cosΛ::Float64`: Cosine of sweep angle
    
    **Outputs:**
    - `t_cap::Float64`: Thickness of cap
    - `Ab_cap::Float64`: Cross-sectional area of cap
"""
function size_cap(σmax, moment, h̄, w̄, h_rms, c, cosΛ)
    t_cap = calc_cap_thickness(σmax, moment, h̄, w̄, h_rms, c, cosΛ)
    Ab_cap = 2 * t_cap * w̄
    return t_cap, Ab_cap
end  # function size_cap

"""
    size_web(τmax, shear, c_perp, web_height)

Calculates area and thickness of wing webs based on maximum shear stress and geometric parameters.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `τmax::Float64`: Maximum allowable shear stress
    - `shear::Float64`: Shear load
    - `c_perp::Float64`: Perpendicular chord length
    - `web_height::Float64`: Height of web
    
    **Outputs:**
    - `t_web::Float64`: Thickness of web
    - `Ab_web::Float64`: Cross-sectional area of web
"""
function size_web(τmax, shear, c_perp, web_height)
    t_web = calc_web_thickness(τmax, shear, c_perp, web_height)
    Ab_web = 2 * t_web * web_height
    return t_web, Ab_web
end  # function size_web

"""
    calc_web_thickness(τmax, shear, c_perp, web_height)

Calculates the required web thickness based on maximum shear stress and loading conditions.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `τmax::Float64`: Maximum allowable shear stress
    - `shear::Float64`: Shear load
    - `c_perp::Float64`: Perpendicular chord length
    - `web_height::Float64`: Height of web
    
    **Outputs:**
    - `t_web::Float64`: Required thickness of web
"""
function calc_web_thickness(τmax, shear, c_perp, web_height)
    t_web = shear / (c_perp^2 * 2 * web_height * τmax)
    return t_web
end  # function calc_web_thickness

"""
    calc_cap_thickness(σmax, moment, h̄, w̄, h_rms, c, cosΛ)

Calculates the required cap thickness based on maximum stress and geometric parameters.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `σmax::Float64`: Maximum allowable stress
    - `moment::Float64`: Bending moment
    - `h̄::Float64`: Height-to-chord ratio
    - `w̄::Float64`: Width-to-chord ratio
    - `h_rms::Float64`: RMS height of sparbox
    - `c::Float64`: Chord length
    - `cosΛ::Float64`: Cosine of sweep angle
    
    **Outputs:**
    - `t_cap::Float64`: Required thickness of cap
"""
function calc_cap_thickness(σmax, moment, h̄, w̄, h_rms, c, cosΛ)
    con = moment * 6h̄ / w̄ * 1 / (c^3 * σmax * cosΛ^4)
    t_cap = 0.5 * (h_rms - ∛(h_rms^3 - con))
    return t_cap
end  # function calc_cap_thickness

"""
    calc_sparbox_internal_area(width, height, t_cap, t_web)

Calculates the internal aera of the sparbox, accounting for the thickness of
the webs and caps.
A = (w - 2tweb)×(h - 2tcap)
"""
function calc_sparbox_internal_area(width, height, t_cap, t_web)
    return (width - 2*t_web)*(height - 2*t_cap)
end  # function calc_internal_area

"""
    add_load!(wing::Wing, load::PointLoad)

Adds a `PointLoad` to the `extra_loads` field of the `Wing` object.
"""
function add_wing_point_load!(wing::Wing, load::PointLoad)  # function added by NILS
    push!(wing.point_loads, load)
end  # function add_wing_point_load!