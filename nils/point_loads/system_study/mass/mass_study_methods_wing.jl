"""
This script calculates the changes in various airframe structural weight components
as a result of applying point loads in different locations and distributions on/across the airframe.
"""

module PointLoadMethods

    export analyse_aircraft, analyse_reference_aircraft

    using Revise
    using TOML
    using StaticArrays
    using Plots
    using TASOPT
    include(TASOPT.__TASOPTindices__)

    function analyse_aircraft(
        ac_segment::String, ac_ref::Any, Weng_single_ref::Any,
        fcs_loc::Any, span_loc::Any, sigma_fcs::Number, wing_frac::Number, nacelle_frac::Number,
    )

        ###

        # Define wing and fuselage point loads
        sigma_fcs_nacelle = sigma_fcs / nacelle_frac

        # Calculate maximum propulsive power throughout mission so far
        T_net_max_ref = ac_ref.pare[ieFe, :] * ac_ref.parg[igneng]  # NILS: see calculate_thrust_from_ROC!() for proof that pare[ieFe] stands for per-engine thrust
        V_0_max_ref = ac_ref.pare[ieu0, :]
        P_prop_max_ref = maximum(T_net_max_ref .* V_0_max_ref)

        # Calculate FCS weight distribution
        if sigma_fcs == -1.0  # GT aircraft
            W_fcs_nacelle_init = Weng_single_ref  # per nacelle
            W_fcs_init = W_fcs_nacelle_init * ac_ref.parg[igneng]  # total on aircraft
        else  # FC aircraft
            W_fcs_nacelle_init = P_prop_max_ref / ac_ref.parg[igneng] / sigma_fcs_nacelle * 9.81  # per nacelle
            W_fcs_init = P_prop_max_ref / sigma_fcs * 9.81  # total on aircraft
        end
        W_fcs_airframe_init = W_fcs_init - ac_ref.parg[igneng] * W_fcs_nacelle_init  # total on airframe
        W_fcs_wing_init = W_fcs_airframe_init * wing_frac / 2  # per wing half
        W_fcs_fuselage_init = W_fcs_airframe_init - 2 * W_fcs_wing_init  # total in fuselage

        ###

        # Write engine point loads to .toml file
        toml_file_path = ""
        toml_ref_file_path = ""
        if ac_segment == "narrowbody"
            global toml_file_path = joinpath(pwd(), "example", "a220100", "a220100_lh2_input.toml")
            # global toml_file_path = joinpath(pwd(), "example", "cryo_input.toml")
        elseif ac_segment == "regional"
            global toml_file_path = joinpath(pwd(), "example", "atr72600", "atr72600_lh2_input.toml")
            # error("message")
        end
        toml_data = TOML.parsefile(toml_file_path)

        # toml_data["Propulsion"]["sigma_fcs_nacelle"] = sigma_fcs_nacelle
        # if span_loc isa Number
        #     toml_data["Propulsion"]["span_loc"] = span_loc
        # end
        # if fcs_loc isa Number
        #     toml_data["Propulsion"]["fcs_loc"] = fcs_loc    
        # end
        # toml_data["Propulsion"]["sigma_fcs"] = sigma_fcs
        # toml_data["Propulsion"]["nacelle_frac"] = nacelle_frac
        # toml_data["Propulsion"]["wing_frac"] = wing_frac
        toml_data["Nils"]["sigma_fcs_nacelle"] = sigma_fcs_nacelle
        # if span_loc isa Number
        toml_data["Nils"]["span_loc"] = span_loc
        # end
        # if fcs_loc isa Number
        toml_data["Nils"]["fcs_loc"] = fcs_loc
        # end
        toml_data["Nils"]["sigma_fcs"] = sigma_fcs
        toml_data["Nils"]["nacelle_frac"] = nacelle_frac
        toml_data["Nils"]["wing_frac"] = wing_frac
        
        println(sigma_fcs_nacelle, " " , sigma_fcs, " ", wing_frac, " ", nacelle_frac, " ", fcs_loc, " ", span_loc)

        open(toml_file_path, "w") do file
            TOML.print(file, toml_data)
        end

        # Read model inputs
        global ac = read_aircraft_model(toml_file_path)  # global needed so updates available in outer scope

        # Add fuselage point load
        ac.fuselage.point_loads = TASOPT.structures.PointLoad[TASOPT.structures.PointLoad()]
        Fz_point_fus_init = -W_fcs_fuselage_init
        if fcs_loc isa String
            _fcs_loc = fcs_loc
        elseif fcs_loc isa Number
            _fcs_loc = Dict("frac_len" => fcs_loc)
        end
        fus_load_init = TASOPT.structures.PointLoad(
            force = SVector(0.0, 0.0, Fz_point_fus_init),  # NILS: previously, I mistakenly multiplied by number of engines here (20.10.2025)
            r = SVector(_fcs_loc, 0.0, 0.0),
            frame = TASOPT.WORLD
        )
        TASOPT.structures.add_fus_point_load!(ac.fuselage, fus_load_init)

        # Add wing point load
        ac.wing.point_loads = TASOPT.structures.PointLoad[TASOPT.structures.PointLoad()]
        Fz_point_wing_init = -W_fcs_wing_init
        if span_loc isa String
            _span_loc = span_loc
        elseif span_loc isa Number
            _span_loc = Dict("frac_span" => span_loc)
        end
        wing_load_init = TASOPT.structures.PointLoad(
            force = SVector(0.0, 0.0, Fz_point_wing_init),
            r = SVector(0.0, _span_loc, 0.0),
            frame = TASOPT.WORLD
        )
        TASOPT.structures.add_wing_point_load!(ac.wing, wing_load_init)

        # Add nacelle point load
        Fz_point_nacelle_init = Weng_single_ref - W_fcs_nacelle_init
        ac.engine.point_load = Fz_point_nacelle_init  # NILS: store point load in engine model

        # Size aircraft
        size_aircraft!(
            ac, iter=100, printiter=false,
            wrlx1=0.5, wrlx2=0.9, wrlx3=0.5,
            # wrlx1=0.5, wrlx2=0.1, wrlx3=0.5,
        )  # disabled default iteration print-outs to console

        return ac, ac_ref

    end

    ###
    
    # function analyse_reference_aircraft(ac_segment::String, N_eng::Int = 2)
    function analyse_reference_aircraft(ac_segment::String)

        # Baseline LH2 aircraft

        # Set all point loads applied via .toml file to 0
        toml_file_path = ""
        toml_ref_file_path = ""
        if ac_segment == "narrowbody"
            global toml_ref_file_path = joinpath(pwd(), "example", "a220100", "a220100_lh2_input_org.toml")
            # global toml_ref_file_path = joinpath(pwd(), "example", "cryo_input_org.toml")
        elseif ac_segment == "regional"
            global toml_ref_file_path = joinpath(pwd(), "example",  "atr72600", "atr72600_lh2_input_org.toml")
            # error("message")
        end
        toml_data = TOML.parsefile(toml_ref_file_path)

        # Read model inputs
        global ac_ref = read_aircraft_model(toml_ref_file_path)

        # # Set point loads applied via PointLoad mutable struct to 0
        # fus_load = TASOPT.structures.PointLoad(
        #     force = SVector(0.0, 0.0, 0.0),
        #     r = SVector(0.0, 0.0, 0.0),
        #     frame = TASOPT.WORLD
        # )
        # TASOPT.structures.add_fus_point_load!(ac_ref.fuselage, fus_load)

        # wing_load = TASOPT.structures.PointLoad(
        #     force = SVector(0.0, 0.0, 0.0),
        #     r = SVector(0.0, 0.0, 0.0),
        #     frame = TASOPT.WORLD
        # )
        # TASOPT.structures.add_wing_point_load!(ac_ref.wing, wing_load)

        # Size aircraft
        size_aircraft!(
            ac_ref, iter=100, printiter=false,
            wrlx1=0.5, wrlx2=0.9, wrlx3=0.5,
        )  # disabled default iteration print-outs to console

        # Print summary data
        # summary(ac_ref)

        # Extract SINGLE engine weight
        Weng_single_ref = ac_ref.parg[igWeng] / ac_ref.parg[igneng]

        return ac_ref, Weng_single_ref

    end
    
end  # module PointLoadMethods
