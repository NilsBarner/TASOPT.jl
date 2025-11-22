using StaticArrays  # added by NILS
using ..structures: PointLoad, add_fus_point_load!, add_wing_point_load!, WORLD  # added by NILS (do NOT use `using TASOPT` inside package files!)
"""
    tfweightwrap!(ac)

General function to estimate and store the weight of a turbofan engine. 
This function is basically a wrapper on tfweight, going from the
basic aircraft inputs to those required by the function and storing the outputs.
      
!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `ac::aircraft`: aircraft object

    **Output:**
    No direct outputs. The `ac` object gets modified with the engine weights.
"""
function tfweightwrap!(ac; Vfcs_is_input::Bool = false)
    """ NILS
    Point loads are added TO engine weight!
    """
    parg = ac.parg
    wing = ac.wing
    neng = parg[igneng]
    
    Weng, Wnace, Webare, W_HXs, Snace1 = tfweight(ac)
    Weng1 = Weng / neng  # NILS: per engine (see `Weng1 = parg[igWeng] / parg[igneng]` in size_aircraft.jl)

    # NILS: Extract custom inputs
    sigma_fcs_nacelle = ac.nils.sigma_fcs_nacelle
    sigma_fcs = ac.nils.sigma_fcs
    wing_frac = ac.nils.wing_frac
    fcs_loc = ac.nils.fcs_loc
    span_loc = ac.nils.span_loc

    # NILS: Calculate maximum propulsive power throughout mission so far
    T_net_history = ac.pare[ieFe, :] * neng  # NILS: see calculate_thrust_from_ROC!() for proof that pare[ieFe] stands for per-engine thrust
    V_0_history = ac.pare[ieu0, :]
    P_prop_max = maximum(T_net_history .* V_0_history)

    # NILS: Calculate FCS weight distribution
    if sigma_fcs == -1.0  # GT aircraft
        W_fcs_nacelle = Weng1  # per nacelle
        W_fcs = W_fcs_nacelle * neng  # total on aircraft
    else  # FC aircraft
        W_fcs_nacelle = P_prop_max / neng / sigma_fcs_nacelle * 9.81  # per nacelle
        W_fcs = P_prop_max / sigma_fcs * 9.81  # total on aircraft
    end
    W_fcs_airframe = W_fcs - neng * W_fcs_nacelle  # total on airframe
    W_fcs_wing = W_fcs_airframe * wing_frac / 2  # per wing half
    W_fcs_fuselage = W_fcs_airframe - 2 * W_fcs_wing  # total in fuselage

    # Add fuselage point load
    ac.fuselage.point_loads = PointLoad[PointLoad()]  # ABSOLUTELY NEEDED, OTHERWISE GETS ADDED WITH EVERY WEIGHT IERATION IN size_aircraft.jl
    Fz_point_fus = -W_fcs_fuselage  # total in fuselage
    if fcs_loc isa Number
        if fcs_loc > 0.0 && fcs_loc < 1.0
            _fcs_loc = Dict("frac_len" => fcs_loc)
        else
            _fcs_loc = fcs_loc
        end
    elseif fcs_loc isa String
        _fcs_loc = fcs_loc
    end
    fus_load = PointLoad(
        force = SVector(0.0, 0.0, Fz_point_fus),
        r = SVector(_fcs_loc, 0.0, 0.0),
        frame = WORLD
    )
    add_fus_point_load!(ac.fuselage, fus_load)

    # Add wing point load
    wing.point_loads = PointLoad[PointLoad()]  # ABSOLUTELY NEEDED, OTHERWISE GETS ADDED WITH EVERY WEIGHT IERATION IN size_aircraft.jl
    Fz_point_wing = -W_fcs_wing  # per wing half
    if span_loc isa Number
        if span_loc > 0.0 && span_loc < 1.0
            _span_loc = Dict("frac_span" => span_loc)
        else
            _span_loc = span_loc
        end
    elseif span_loc isa String
        _span_loc = span_loc
    end
    wing_load = PointLoad(
        force = SVector(0.0, 0.0, Fz_point_wing),
        r = SVector(0.0, _span_loc, 0.0),
        frame = WORLD
    )
    add_wing_point_load!(ac.wing, wing_load)
    
    # Add nacelle point load
    Fz_point_nacelle = Weng1 - W_fcs_nacelle  # per nacelle
    ac.engine.point_load = Fz_point_nacelle  # store point load in engine model
    parg[igWeng] = Weng - neng * Fz_point_nacelle  # second term added by NILS
    parg[igWebare] = Webare
    parg[igWnace] = Wnace
    parg[igWHXs] = W_HXs

    # set new nacelle area / reference area  fraction fSnace
    S = wing.layout.S

    if Vfcs_is_input == true
        # NILS: correct nacelle volume (and thus surface area) by FCS volume
        dfan = ac.parg[igdfan]
        rSnace = parg[igrSnace]
        Vfcs = ac.nils.V_fcs_nacelle / neng  # parg[igVfcsnac]
        lfcs = Vfcs / (pi * dfan^2)  # assume available FCS frontal area to equal fan hub area
        Snace = Snace1 * (rSnace + lfcs) / rSnace * neng  # second term in brackets added by NILS
    elseif Vfcs_is_input == false
        Snace = Snace1 * neng
        lfcs = 0.0  # do not add length increment
    end

    # NILS: the original tfweightwrap.jl filed contained a few repeated lines here which I removed
    fSnace = Snace / S
    parg[igfSnace] = fSnace
    lnace = parg[igdfan] * parg[igrSnace] * 0.15 + lfcs  # last term added by NILS
    parg[iglnace] = lnace

    if Vfcs_is_input == true
        nothing
    elseif Vfcs_is_input == false
        # NILS: calculate and log total available nacelle volume for FCS
        # assuming FCS volume to be a truncated cone. The below expression
        # for the front and aft radius was taken from stickfig() in output_plots.jl
        # see lines following `# Aft fan outline` and is more qualitative than accurate.
        dfan = ac.parg[igdfan]
        HTR_f = parg[igHTRf]
        V_fcs_av_nace = 1/3 * pi * lnace * ((dfan/2 * HTR_f)^2 + (dfan/4 * HTR_f)^2 + (dfan/2 * HTR_f) * (dfan/4 * HTR_f))
        # println("lnace, dfan, HTR_f, V_fcs_av_nace = ", lnace, ", ", dfan, ", ", HTR_f, ", ", V_fcs_av_nace)
        # parg[igVfcsavnacetot] = neng * V_fcs_av_nace  # ALL engines
        ac.nils.V_fcs_nacelle = neng * V_fcs_av_nace  # ALL engines
    end

end