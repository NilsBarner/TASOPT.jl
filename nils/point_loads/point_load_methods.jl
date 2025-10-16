"""
This script calculates the changes in various airframe structural weight components
as a resuls of applying point loads in different locations on the airframe.
It is being called from different execution scripts to analyse various load cases.
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
        ac_segment::String, case_idx::Int, _ac_ref::Any, _Weng_single_ref::Any, indep_var::Any, fcs_loc::Any, N_eng::Int,
        sigma_gt::Number, sigma_fc_wing_min::Number, sigma_fc_wing_max::Number, sigma_fc_min::Number, sigma_fc_max::Number,
    )

        # Size reference aircraft if not provided (only needed for case 5)
        if _ac_ref === nothing
            ac_ref, Weng_single_ref = analyse_reference_aircraft(ac_segment, N_eng)
        else
            ac_ref, Weng_single_ref = _ac_ref, _Weng_single_ref
        end

        # Calculate FCS weight limits
        W_fc_wing_min = Weng_single_ref * sigma_gt / sigma_fc_wing_max
        W_fc_wing_max = Weng_single_ref * sigma_gt / sigma_fc_wing_min
        W_fc_min = Weng_single_ref * sigma_gt / sigma_fc_max
        W_fc_max = Weng_single_ref * sigma_gt / sigma_fc_min

        # Set point loads by case
        if in(case_idx, (1, 61))
            W_fc = indep_var
            W_fc_wing = W_fc
            global Fz_point_fus = 0.0
            global Fz_point_wing = -(W_fc_wing - Weng_single_ref)  # -ve
        elseif in(case_idx, (21, 22, 621, 622))
            W_fc = indep_var
            W_fc_wing = W_fc_wing_min
            W_fc_fus = W_fc - W_fc_wing_min
            global Fz_point_wing = -(W_fc_wing - Weng_single_ref)  # +ve
            global Fz_point_fus = -W_fc_fus  # -ve
        elseif in(case_idx, (31, 32, 631, 632))
            W_fc_wing = indep_var
            W_fc_fus = W_fc_wing_max - W_fc_wing
            global Fz_point_fus = -W_fc_fus  # -ve
            global Fz_point_wing = -(W_fc_wing - Weng_single_ref)  # -ve
        elseif in(case_idx, (4, 64))
            # fcs_loc = indep_var  # DO NOT USE as fcs_loc must be of type Dict for the logic in fuseW.jl to work
            W_fc = W_fc_min
            W_fc_wing = W_fc_wing_min
            W_fc_fus = W_fc - W_fc_wing_min
            global Fz_point_wing = -(W_fc_wing - Weng_single_ref)  # +ve
            global Fz_point_fus = -W_fc_fus  # -ve
        elseif in(case_idx, (5, 65))
            N_eng = indep_var
            W_fc = W_fc_min
            W_fc_wing = W_fc
            global Fz_point_fus = 0.0
            global Fz_point_wing = -(W_fc_wing - Weng_single_ref)  # -ve
        end

        # Write engine point loads to .toml file
        if ac_segment == "narrowbody"
            toml_data = TOML.parsefile(joinpath(pwd(), "example", "cryo_input.toml"))
        elseif ac_segment == "regional"
            toml_data = TOML.parsefile(joinpath(pwd(), "nils", "point_loads", "default_regional_cryo.toml"))
        end
        toml_data["Propulsion"]["number_of_engines"] = N_eng
        toml_data["Propulsion"]["Weight"]["custom_weight_delta"] = -Fz_point_wing  # negative sign accounts for positive weights
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

        # Add fuselage point loads to mutable structure
        fus_load = TASOPT.structures.PointLoad(
            force = SVector(0.0, 0.0, toml_data["Propulsion"]["number_of_engines"] * Fz_point_fus),
            r = SVector(fcs_loc, 0.0, 0.0),
            frame = TASOPT.WORLD
        )
        TASOPT.structures.add_fus_point_load!(ac.fuselage, fus_load)

        # Size aircraft
        size_aircraft!(ac, iter=50, printiter=false)  # disabled default iteration print-outs to console

        return ac, ac_ref

    end

    ###

    function analyse_reference_aircraft(ac_segment::String, N_eng::Int = 2)

        # Baseline LH2 aircraft

        # Set all point loads applied via .toml file to 0
        if ac_segment == "narrowbody"
            toml_data = TOML.parsefile(joinpath(pwd(), "example", "cryo_input.toml"))
        elseif ac_segment == "regional"
            toml_data = TOML.parsefile(joinpath(pwd(), "nils", "point_loads", "default_regional_cryo.toml"))
        end
        toml_data["Propulsion"]["number_of_engines"] = N_eng
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

        return ac_ref, Weng_single_ref

    end

end  # module PointLoadMethods
