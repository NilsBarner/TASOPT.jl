using Revise
using TOML
using StaticArrays
using Plots
using TASOPT
include(TASOPT.__TASOPTindices__)

include(joinpath(pwd(), "nils/point_loads/point_load_plots.jl"))
using .PointLoadPlots

using XLSX
using DataFrames
using Tables

using CSV
include(joinpath(pwd(), "nils/point_loads/system_study/volume/wingbox_coordinate_compiler.jl"))

###

# Top-level settings
ac_segment = "narrowbody"  # "narrowbody" or "regional"
##### NILS: RESET "REGIONAL" - CURRENTLY DOES NOT RUN, PRESUMABLY BECAUSE OF INFEASIBLE INPUT FROM MASS-EXPERIMENTS #####
N_points = 19

###

function analyse_reference_aircraft(ac_segment::String)

    # Baseline LH2 aircraft

    # Set all point loads applied via .toml file to 0
    if ac_segment == "narrowbody"
        toml_data = TOML.parsefile(joinpath(pwd(), "example", "cryo_input.toml"))
    elseif ac_segment == "regional"
        toml_data = TOML.parsefile(joinpath(pwd(), "nils", "point_loads", "default_regional_cryo.toml"))
    end
    toml_data["Propulsion"]["number_of_engines"] = 2
    toml_data["Propulsion"]["Turbomachinery"]["HTR_fan"] = 0.3
    toml_data["Propulsion"]["V_fcs_nacelle"] = 0.0
    toml_data["Wing"]["V_fcs_wing"] = 0.0
    toml_data["Fuselage"]["V_fcs_fuselage"] = 0.0
    toml_data["Fuselage"]["l_fcs_fuselage"] = 0.0
    toml_data["Fuselage"]["Geometry"]["unit_load_device"] = "LD3-45"  # NILS: "LD3" is too big (see src/utils/constants.jl)
    toml_data["Fuselage"]["Geometry"]["theta_floor"] = 0.0  # [deg]
    toml_data["Fuselage"]["fcs_fuselage_location"] = 0.0  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl"
    toml_data["Wing"]["AR"] = 10.1
    toml_data["Wing"]["tdivc_scale"] = 1.0
    toml_data["Wing"]["TR_scale"] = 1.0
    toml_data["Wing"]["has_strut"] = false
    toml_data["Fuselage"]["Geometry"]["radius"] = 100 * 0.0254  # [m]
    toml_data["Propulsion"]["Weight"]["custom_weight_delta"] = 0.0  # point load at engine location on wing
    if ac_segment == "narrowbody"
        open(joinpath(pwd(), "example", "cryo_input.toml"), "w") do file
            TOML.print(file, toml_data)
        end
    elseif ac_segment == "regional"
        open(joinpath(pwd(), "nils", "point_loads", "default_regional_cryo.toml"), "w") do file
            TOML.print(file, toml_data)
        end
    end

    # Read model inputs
    if ac_segment == "narrowbody"
        global ac_ref = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))
    elseif ac_segment == "regional"
        global ac_ref = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../nils/point_loads/default_regional_cryo.toml"))
    end

    # Set point loads applied via PointLoad mutable struct to 0
    fus_load = TASOPT.structures.PointLoad(
        force = SVector(0.0, 0.0, 0.0),
        r = SVector(0.0, 0.0, 0.0),
        frame = TASOPT.WORLD
    )
    TASOPT.structures.add_fus_point_load!(ac_ref.fuselage, fus_load)

    # Size aircraft
    size_aircraft!(ac_ref, iter=50, printiter=false)  # disabled default iteration print-outs to console

    # Print summary data
    # summary(ac_ref)

    # Extract SINGLE engine weight
    Weng_single_ref = ac_ref.parg[igWeng] / ac_ref.parg[igneng]

    ###

    # Calculate wing box centroid location
    xw, yw, zw = compile_wingbox_coordinates(ac_ref)
    x_centroid_wing = sum(xw) / length(xw)
    y_centroid_wing = sum(yw) / length(yw)
    z_centroid_wing = sum(zw) / length(zw)

    # Calculate engine centroid location (adopted from around line 933 in src/sizing/wsize.jl)
    b = ac_ref.wing.layout.span
    bo = ac_ref.wing.layout.root_span
    ηs = ac_ref.wing.layout.ηs
    dy = 2 * ac_ref.parg[igdfan] # space to leave near wing root and tip [m]
    if toml_data["Propulsion"]["number_of_engines"] == 2
        yi = ηs * b / 2
    else
        yi = LinRange(bo / 2 + dy, b / 2 * 3 / 4, Int(toml_data["Propulsion"]["number_of_engines"] / 2))
    end
    y_centroid_nacelles = sum(yi) / length(yi)

    df = DataFrame(
        index = 0.0,
        fcs_loc = "NA",
        radius = toml_data["Fuselage"]["Geometry"]["radius"],
        AR = toml_data["Wing"]["AR"],
        tdivc_scale = toml_data["Wing"]["tdivc_scale"],
        N_eng = toml_data["Propulsion"]["number_of_engines"],
        HTR_f = toml_data["Propulsion"]["Turbomachinery"]["HTR_fan"],
        Vspec = 1e30,
        fcs_fuselage_location = 0.0,
        has_strut = toml_data["Wing"]["has_strut"],
        CDS = ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1],
        WMTO = ac_ref.parg[igWMTO],
        Vol_wing = 2*(ac_ref.wing.center.volume + ac_ref.wing.inboard.volume + ac_ref.wing.outboard.volume),
        PFEI = ac_ref.parm[imPFEI],
        seats_abreast = ac_ref.fuselage.cabin.seats_abreast_main,
        Vol_nacelle = ac_ref.parg[igVfcsavnacetot],
        y_centroid_wing = y_centroid_wing,
        L_fuse = ac_ref.fuselage.layout.l_cabin_cylinder,
        y_centroid_nacelles = y_centroid_nacelles,
        span = b,
    )

    # println("CDS =", ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1])
    # println("WMTO =", ac_ref.parg[igWMTO])

    # CSV.write("volume_results_narrowbody_ref_newest.csv", df)  # fasts

    return ac_ref, Weng_single_ref

end

###

function run_study()

    # Model reference aircraft
    ac_ref, Weng_single_ref = analyse_reference_aircraft(ac_segment)
    error("message")

    ###


    P_fcs = 20e6  # [W] total FCS power

    # Set all point loads applied via .toml file to 0
    if ac_segment == "narrowbody"
        toml_data = TOML.parsefile(joinpath(pwd(), "example", "cryo_input.toml"))
    elseif ac_segment == "regional"
        toml_data = TOML.parsefile(joinpath(pwd(), "nils", "point_loads", "default_regional_cryo.toml"))
    end
    radius_ref = toml_data["Fuselage"]["Geometry"]["radius"]
    AR_ref = toml_data["Wing"]["AR"]
    tdivc_scale_ref = toml_data["Wing"]["tdivc_scale"]
    N_eng_ref = toml_data["Propulsion"]["number_of_engines"]
    HTR_f_ref = toml_data["Propulsion"]["Turbomachinery"]["HTR_fan"]
    Vspec_ref = 1e30
    fcs_fuselage_location_ref = toml_data["Fuselage"]["fcs_fuselage_location"]
    has_strut_ref = false

    fcs_locs = ["nacelle", "wing", "fuselage"]
    # fcs_locs = ["wing"]
    # fcs_locs = ["nacelle"]
    # fcs_locs = ["fuselage"]
    radius_range = collect(range(80, 130, length=N_points)) * 0.0254  # [m]
    AR_range = collect(range(6, 18, length=N_points))
    tdivc_scale_range = collect(range(0.75, 1.25, length=N_points))
    Neng_range = [2, 4, 6, 8]
    HTRf_range = collect(range(0.2, 0.85, length=N_points))
    # Vspec_range = range(3, 5.5, length=N_points) * 1e3 * 1e3  # [m^3]  Table 2 in FZO-PPN-REP-0031 - Fuel Cells Technical Report (STACK)
    # Vspec_range = range(0.2, 0.4, length=N_points) * 1e3 * 1e3  # [m^3]  Table 2 in FZO-PPN-REP-0031 - Fuel Cells Technical Report (SYSTEM)
    Vspec_range = collect(range(0.2, 1.0, length=N_points)) * 1e3 * 1e3  # NILS' estimate
    fcs_fuselage_locations = [0.0, 1.0]
    strut_bools = [false, true]

    # --- insert current single values while keeping ordering and uniqueness ---
    radius_range = sort(unique(vcat(Float64.(radius_range), Float64(radius_ref))))
    AR_range = sort(unique(vcat(Float64.(AR_range), Float64(AR_ref))))
    tdivc_scale_range = sort(unique(vcat(Float64.(tdivc_scale_range), Float64(tdivc_scale_ref))))
    HTRf_range = sort(unique(vcat(Float64.(HTRf_range), Float64(HTR_f_ref))))
    Vspec_range = sort(unique(vcat(Float64.(Vspec_range), Float64(Vspec_ref))))
    # fcs_locs = ["nacelle"]
    # radius_range = [radius_ref]
    # AR_range = [AR_ref]
    # tdivc_scale_range = [tdivc_scale_ref]
    # Neng_range = [2]
    # HTRf_range = [HTR_f_ref]
    # Vspec_range = [Vspec_ref]
    # fcs_fuselage_locations = [0.0]

    # result vectors
    index_vec = Vector{NTuple{9,Int}}()
    fcs_loc_vec = String[]
    radius_vec = Float64[]
    AR_vec = Float64[]
    tdivc_vec = Float64[]
    Neng_vec = Int[]
    HTRf_vec = Float64[]
    Vspec_vec = Float64[]
    fcs_fuselage_loc_vec = Float64[]
    strut_bool_vec = Bool[]

    CDS_vec = Float64[]
    WMTO_vec = Float64[]
    Vol_wing_vec = Float64[]
    PFEI_vec = Float64[]
    seats_abreast_vec = Int[]
    Vol_nacelle_vec = Float64[]
    y_centroid_wing_vec = Float64[]
    L_fuse_vec = Float64[]
    y_centroid_nacelles_vec = Float64[]
    span_vec = Float64[]

    for (h, fcs_loc) in enumerate(fcs_locs)
        for (i, radius) in enumerate(radius_range)
            for (j, AR) in enumerate(AR_range)
                for (k, tdivc_scale) in enumerate(tdivc_scale_range)
                    for (l, N_eng) in enumerate(Neng_range)
                        for (m, HTR_f) in enumerate(HTRf_range)
                            for (n, Vspec) in enumerate(Vspec_range)
                                for (o, fcs_fuselage_location) in enumerate(fcs_fuselage_locations)
                                    for (p, has_strut) in enumerate(strut_bools)

                                        # your condition that decides whether to compute
                                        if (h == 1 && i == 1 && j == 1 && k == 1 && n == 1 && o == 1 && p == 1) ||
                                           (h == 2 && i == 1 && l == 1 && m == 1 && n == 1 && o == 1) ||
                                           (h == 3 && j == 1 && k == 1 && l == 1 && m == 1 && p == 1)

                                        # if (i == 1 && l == 1 && m == 1 && n == 1 && o == 1)  # wing
                                        # if (i == 1 && j == 1 && k == 1 && n == 1 && o == 1)  # nacelle
                                        # if (j == 1 && k == 1 && l == 1 && m == 1)  # fuselage

                                            if (h == 1 && i == 1 && j == 1 && k == 1 && n == 1 && o == 1 && p == 1)
                                                # radius = 100 * 0.0254  # [m]
                                                # AR = 10.1
                                                # tdivc_scale = 1.0
                                                # Vspec = 0.5e6
                                                # fcs_fuselage_location = 0.0
                                                # has_strut = false
                                                _h = h
                                                _i, radius = findfirst(==(radius_ref), radius_range), radius_ref
                                                _j, AR = findfirst(==(AR_ref), AR_range), AR_ref
                                                _k, tdivc_scale = findfirst(==(tdivc_scale_ref), tdivc_scale_range), tdivc_scale_ref
                                                _l = l
                                                _m = m
                                                _n, Vspec = findfirst(==(Vspec_ref), Vspec_range), Vspec_ref
                                                _o, fcs_fuselage_location = findfirst(==(fcs_fuselage_location_ref), fcs_fuselage_locations), fcs_fuselage_location_ref
                                                _p, has_strut = findfirst(==(has_strut_ref), strut_bools), has_strut_ref
                                                # fixed println: supply commas between arguments
                                                # println("n = ", n, " Vspec = ", Vspec)
                                                # error("message")

                                            elseif (h == 2 && i == 1 && l == 1 && m == 1 && n == 1 && o == 1)
                                                # radius = 100 * 0.0254  # [m]
                                                # N_eng = 2
                                                # HTR_f = 0.3
                                                # Vspec = 0.5e6
                                                # fcs_fuselage_location = 0.0
                                                _h = h
                                                _i, radius = findfirst(==(radius_ref), radius_range), radius_ref
                                                _j = j
                                                _k = k
                                                _l, N_eng = findfirst(==(N_eng_ref), Neng_range), N_eng_ref
                                                _m, HTR_f = findfirst(==(HTR_f_ref), HTRf_range), HTR_f_ref
                                                _n, Vspec = findfirst(==(Vspec_ref), Vspec_range), Vspec_ref
                                                _o, fcs_fuselage_location = findfirst(==(fcs_fuselage_location_ref), fcs_fuselage_locations), fcs_fuselage_location_ref
                                                _p = p
                                            elseif (h == 3 && j == 1 && k == 1 && l == 1 && m == 1 && p == 1)
                                                # AR = 10.1
                                                # tdivc_scale = 1.0
                                                # N_eng = 2
                                                # HTR_f = 0.3
                                                # has_strut = false
                                                _h = h
                                                _i = i
                                                _j, AR = findfirst(==(AR_ref), AR_range), AR_ref
                                                _k, tdivc_scale = findfirst(==(tdivc_scale_ref), tdivc_scale_range), tdivc_scale_ref
                                                _l, N_eng = findfirst(==(N_eng_ref), Neng_range), N_eng_ref
                                                _m, HTR_f = findfirst(==(HTR_f_ref), HTRf_range), HTR_f_ref
                                                _n = n
                                                _o = o
                                                _p, has_strut = findfirst(==(has_strut_ref), strut_bools), has_strut_ref
                                            end

                                            V_fcs = P_fcs / Vspec  # [m^3] total FCS volume

                                            ###

                                            # Write engine point loads to .toml file
                                            if ac_segment == "narrowbody"
                                                toml_data = TOML.parsefile(joinpath(pwd(), "example", "cryo_input.toml"))
                                            elseif ac_segment == "regional"
                                                toml_data = TOML.parsefile(joinpath(pwd(), "nils", "point_loads", "default_regional_cryo.toml"))
                                            end
                                            toml_data["Propulsion"]["number_of_engines"] = N_eng
                                            toml_data["Propulsion"]["Turbomachinery"]["HTR_fan"] = HTR_f
                                            if fcs_loc == "nacelle"
                                                toml_data["Propulsion"]["V_fcs_nacelle"] = V_fcs / N_eng
                                                toml_data["Wing"]["V_fcs_wing"] = 0.0
                                                toml_data["Fuselage"]["V_fcs_fuselage"] = 0.0
                                                toml_data["Fuselage"]["l_fcs_fuselage"] = 0.0
                                                toml_data["Wing"]["has_strut"] = false
                                            elseif fcs_loc == "wing"                                                
                                                toml_data["Propulsion"]["V_fcs_nacelle"] = 0.0
                                                toml_data["Wing"]["V_fcs_wing"] = V_fcs
                                                toml_data["Fuselage"]["V_fcs_fuselage"] = 0.0
                                                toml_data["Fuselage"]["l_fcs_fuselage"] = 0.0
                                                toml_data["Wing"]["has_strut"] = has_strut
                                            elseif fcs_loc == "fuselage"
                                                toml_data["Propulsion"]["V_fcs_nacelle"] = 0.0
                                                toml_data["Wing"]["V_fcs_wing"] = 0.0
                                                # toml_data["Fuselage"]["V_fcs_fuselage"] = V_fcs
                                                toml_data["Fuselage"]["l_fcs_fuselage"] = l_fcs
                                                toml_data["Fuselage"]["fcs_fuselage_location"] = fcs_fuselage_location  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl"
                                                toml_data["Wing"]["has_strut"] = false
                                            end
                                            toml_data["Wing"]["AR"] = AR
                                            toml_data["Wing"]["tdivc_scale"] = tdivc_scale
                                            toml_data["Wing"]["TR_scale"] = 1.0
                                            toml_data["Fuselage"]["Geometry"]["radius"] = radius
                                            if ac_segment == "narrowbody"
                                                open(joinpath(pwd(), "example", "cryo_input.toml"), "w") do file
                                                    TOML.print(file, toml_data)
                                                end
                                            elseif ac_segment == "regional"
                                                open(joinpath(pwd(), "nils", "point_loads", "default_regional_cryo.toml"), "w") do file
                                                    TOML.print(file, toml_data)
                                                end
                                            end

                                            # Read model inputs
                                            if ac_segment == "narrowbody"
                                                global ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../example/cryo_input.toml"))  # global needed so updates available in outer scope
                                            elseif ac_segment == "regional"
                                                global ac = read_aircraft_model(joinpath(TASOPT.__TASOPTroot__, "../nils/point_loads/default_regional_cryo.toml"))  # global needed so updates available in outer scope
                                            end 
                                        
                                            # (wrap in try/catch as before)
                                            try
                                                # size_aircraft! and get ac
                                                size_aircraft!(ac, iter=50, printiter=false)

                                                # Calculate wing box centroid location
                                                xw, yw, zw = compile_wingbox_coordinates(ac)
                                                x_centroid_wing = sum(xw) / length(xw)
                                                y_centroid_wing = sum(yw) / length(yw)
                                                z_centroid_wing = sum(zw) / length(zw)

                                                # Calculate engine centroid location (adopted from around line 933 in src/sizing/wsize.jl)
                                                b = ac.wing.layout.span
                                                bo = ac.wing.layout.root_span
                                                ηs = ac.wing.layout.ηs
                                                dy = 2 * ac.parg[igdfan] # space to leave near wing root and tip [m]
                                                if N_eng == 2
                                                    yi = ηs * b / 2
                                                else
                                                    yi = LinRange(bo / 2 + dy, b / 2 * 3 / 4, Int(N_eng / 2))
                                                end
                                                y_centroid_nacelles = sum(yi) / length(yi)

                                                # save results immediately
                                                push!(index_vec, (_h,_i,_j,_k,_l,_m,_n,_o,_p))
                                                push!(fcs_loc_vec, fcs_loc)
                                                push!(radius_vec, radius)
                                                push!(AR_vec, AR)
                                                push!(tdivc_vec, tdivc_scale)
                                                push!(Neng_vec, N_eng)
                                                push!(HTRf_vec, HTR_f)
                                                push!(Vspec_vec, Vspec)
                                                push!(fcs_fuselage_loc_vec, fcs_fuselage_location)
                                                push!(strut_bool_vec, has_strut)

                                                push!(CDS_vec, ac.wing.layout.S * ac.para[iaCD, ipcruise1])
                                                push!(WMTO_vec, ac.parg[igWMTO])
                                                push!(Vol_wing_vec, 2*(ac.wing.center.volume + ac.wing.inboard.volume + ac.wing.outboard.volume))
                                                push!(PFEI_vec, ac.parm[imPFEI])
                                                push!(seats_abreast_vec, ac.fuselage.cabin.seats_abreast_main)
                                                push!(Vol_nacelle_vec, ac.parg[igVfcsavnacetot])
                                                push!(y_centroid_wing_vec, y_centroid_wing)
                                                push!(L_fuse_vec, ac.fuselage.layout.l_cabin_cylinder)
                                                push!(y_centroid_nacelles_vec, y_centroid_nacelles)
                                                push!(span_vec, b)

                                                # println("CDS =", ac.wing.layout.S * ac.para[iaCD, ipcruise1])
                                                # println("WMTO =", ac.parg[igWMTO])
                                                # error("message")

                                            catch e
                                                # @warn "failed case" fcs_loc=fcs_loc radius=radius AR=AR exception=e
                                                println("fcs_loc: $fcs_loc, radius: $radius, AR: $AR, tdivc_scale: $tdivc_scale, N_eng: $N_eng, HTR_f: $HTR_f, Vspec: $Vspec, fcs_fuselage_location: $fcs_fuselage_location", " has_strut: $has_strut")
                                            end

                                        end
                                    
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    df = DataFrame(
      index = index_vec,
      fcs_loc = fcs_loc_vec,
      radius = radius_vec,
      AR = AR_vec,
      tdivc_scale = tdivc_vec,
      N_eng = Neng_vec,
      HTR_f = HTRf_vec,
      Vspec = Vspec_vec,
      fcs_fuselage_location = fcs_fuselage_loc_vec,
      has_strut = strut_bool_vec,
      CDS = CDS_vec,
      WMTO = WMTO_vec,
      Vol_wing = Vol_wing_vec,
      PFEI = PFEI_vec,
      seats_abreast = seats_abreast_vec,
      Vol_nacelle = Vol_nacelle_vec,
      y_centroid_wing = y_centroid_wing_vec,
      L_fuse = L_fuse_vec,
      y_centroid_nacelles = y_centroid_nacelles_vec,
      span = span_vec,
    )

    CSV.write("volume_results_narrowbody_newestttt.csv", df)  # fasts
    # CSV.write("volume_results_narrowbody_wing_newest.csv", df)  # fast
    # CSV.write("volume_results_narrowbody_nacelle_newest.csv", df)  # fast
    # CSV.write("volume_results_narrowbody_fuselage_newest.csv", df)  # fast
    return df
end

df = run_study()

