"""
This script calculates the changes in various airframe structural weight components
as a resuls of applying point loads in different locations on the airframe.
It is being called from different execution scripts to analyse various load cases:

1)  Variable weight (2-4 kW/kg) in fixed location on wing (coincident with 2 engines)
    Variable weight (2-4 kW/kg) in fixed location within fuselage
21) At LH2 tank centre
22) At wing box centre
    Variable split (0-1) of fixed weight (2 kW/kg) between wing and fuselage with fuselage load
31) At LH2 tank centre
32) At wing box centre
4)  Variable location (aft end nose cone to front end tail cone) of fixed weight (2 kW/kg)
5)  Distribution of fixed weight (2 kW/kg) across variable number of engines (2-8)
"""

using Revise
using TOML
using StaticArrays
using Plots
using TASOPT
include(TASOPT.__TASOPTindices__)

include(joinpath(pwd(), "nils/point_loads/point_load_plots.jl"))
using .PointLoadPlots

include(joinpath(pwd(), "nils/point_loads/point_load_methods.jl"))
using .PointLoadMethods

###

ac_segment = "narrowbody"  # "narrowbody" or "regional"
N_points = 4

###

# Model reference aircraft
ac_ref, Weng_single_ref = PointLoadMethods.analyse_reference_aircraft(ac_segment)

###

# Data

sigma_gt = 8e3  # W/kg based on CFM Leap in GT_baseline_slides.pptx
sigma_fc_wing_min = 2e3  # W/kg (entire FCS in nacelle)
sigma_fc_wing_max = 12e3  # W/kg (only motor + inverter in nacelle)
sigma_fc_min = 2e3  # W/kg (minimum for entire FCS in nacelle)
sigma_fc_max = 4e3  # W/kg (maximum for entire FCS in nacelle)

# Calculate FCS weight limits
W_fc_wing_min = Weng_single_ref * sigma_gt / sigma_fc_wing_max
W_fc_wing_max = Weng_single_ref * sigma_gt / sigma_fc_wing_min
W_fc_min = Weng_single_ref * sigma_gt / sigma_fc_max
W_fc_max = Weng_single_ref * sigma_gt / sigma_fc_min

###

# Define results dictionary for storage (pick one line)

# case_idx_list = [1, 21, 22, 31, 32, 4, 5]
# case_idx_list = [1, 21, 22]
# case_idx_list = [31, 32]
# case_idx_list = [4]
# case_idx_list = [5]
# case_idx_list = [1, 61]
case_idx_list = [1, 21, 22, 31, 32, 4, 5, 61, 621, 622, 631, 632, 64, 65]

result_dict = Dict(
    i => Dict(
        "indep_var_list" => [],
        "indep_var_list_1" => [],
        "indep_var_list_2" => [],
        "indep_var_grid" => Array{Any}(undef, N_points, N_points),
        # "indep_var_grid" => fill(NaN, N_points, 4),
        "ac_list" => [],
        "ac_grid" => Array{Any}(undef, N_points, N_points),
        "ac_ref" => Union{Nothing, aircraft},
        "Wpointload_fuselage" => [],
    ) for i in case_idx_list
)

###

# Iterate over cases
for (h, case_idx) in enumerate(case_idx_list)

    global indep_var_list = []
    global indep_var_list_1 = []
    global indep_var_list_2 = []
    # global indep_var_grid = fill(NaN, N_points, N_points)
    global indep_var_grid = fill(NaN, N_points, 4)

    if in(case_idx, (61, 621, 622, 631, 632, 64, 65))
        global ac_segment = "regional"

        # Model reference aircraft
        global ac_ref, Weng_single_ref = PointLoadMethods.analyse_reference_aircraft(ac_segment)
        
        # Calculate FCS weight limits
        global W_fc_wing_min = Weng_single_ref * sigma_gt / sigma_fc_wing_max
        global W_fc_wing_max = Weng_single_ref * sigma_gt / sigma_fc_wing_min
        global W_fc_min = Weng_single_ref * sigma_gt / sigma_fc_max
        global W_fc_max = Weng_single_ref * sigma_gt / sigma_fc_min
    end

    if (case_idx == 1 || case_idx == 61)
        indep_var_list = range(W_fc_min, stop=W_fc_max, length=N_points)  # entire FCS weight varies
        fcs_loc = 0.0  # value not of importance
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
        N_points_1, N_points_2 = N_points, N_points

    elseif (case_idx == 21 || case_idx == 621)
        indep_var_list = range(W_fc_min, stop=W_fc_max, length=N_points)  # entire FCS weight varies
        fcs_loc = "apu"
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
        N_points_1, N_points_2 = N_points, N_points

    elseif (case_idx == 22 || case_idx == 622)
        indep_var_list = range(W_fc_min, stop=W_fc_max, length=N_points)  # entire FCS weight varies
        fcs_loc = "eng"
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
        N_points_1, N_points_2 = N_points, N_points

    elseif (case_idx == 31 || case_idx == 631)

        sigma_fc_range = range(2e3, stop=4e3, length=N_points)
        indep_var_list_1 = sigma_fc_range
        wing_frac_range = range(0.0, stop=1.0, length=N_points)
        indep_var_list_2 = wing_frac_range
        indep_var_grid = zeros(N_points, N_points)
        for (i, sigma_fc) in enumerate(sigma_fc_range)
            for (j, wing_frac) in enumerate(wing_frac_range)
                sigma_fc_wing = sigma_fc_wing_max - wing_frac * (sigma_fc_wing_max - sigma_fc)
                W_fc_wing = Weng_single_ref * sigma_gt / sigma_fc_wing
                indep_var_grid[i, j] = W_fc_wing
            end
        end

        fcs_loc = "apu"
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
        N_points_1, N_points_2 = N_points, N_points

    elseif (case_idx == 32 || case_idx == 632)

        sigma_fc_range = range(2e3, stop=4e3, length=N_points)
        indep_var_list_1 = sigma_fc_range
        wing_frac_range = range(0.0, stop=1.0, length=N_points)
        indep_var_list_2 = wing_frac_range
        indep_var_grid = zeros(N_points, N_points)
        for (i, sigma_fc) in enumerate(sigma_fc_range)
            for (j, wing_frac) in enumerate(wing_frac_range)
                sigma_fc_wing = sigma_fc_wing_max - wing_frac * (sigma_fc_wing_max - sigma_fc)
                W_fc_wing = Weng_single_ref * sigma_gt / sigma_fc_wing
                indep_var_grid[i, j] = W_fc_wing
            end
        end

        fcs_loc = "eng"
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
        N_points_1, N_points_2 = N_points, N_points

    elseif (case_idx == 4 || case_idx == 64)
        indep_var_list_1 = range(sigma_fc_min, stop=sigma_fc_max, length=N_points)
        indep_var_list_2 = range(0, stop=1, length=N_points)  # location of point load varies
        fcs_loc = 0.0  # value not of importance
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
        N_points_1, N_points_2 = N_points, N_points

    elseif (case_idx == 5 || case_idx == 65)
        indep_var_list_1 = range(sigma_fc_min, stop=sigma_fc_max, length=N_points)
        println("indep_var_list_1: ", indep_var_list_1)
        indep_var_list_2 = [2, 4, 6, 8]
        println("indep_var_list_2: ", indep_var_list_2)
        fcs_loc = 0.0  # value not of importance
        _ac_ref, _Weng_single_ref = nothing, nothing
        N_eng = 0  # value not of importance
        N_points_1, N_points_2 = N_points, length(indep_var_list_2)
    end

    # Output lists
    global ac_list = []
    global Wpointload_fuselage_list = Float64[]

    global ac_grid = Array{Any}(undef, N_points_1, N_points_2)
    global Wpointload_fuselage_grid = Array{Float64}(undef, N_points_1, N_points_2)

    if (case_idx == 1 || case_idx == 61 || case_idx == 21 || case_idx == 621 || case_idx == 22 || case_idx == 622)

        sigma_fc = sigma_fc_min

        # Iterative over independent variable
        for (i, indep_var) in enumerate(indep_var_list)

            # try  # FOR DEBUGGING, comment try-catch block

            # Size aircraft
            global ac, ac_ref = PointLoadMethods.analyse_aircraft(
                ac_segment, case_idx, _ac_ref, _Weng_single_ref, indep_var, fcs_loc, N_eng,
                sigma_gt, sigma_fc_wing_min, sigma_fc_wing_max, sigma_fc,  # sigma_fc_min, sigma_fc_max
            )  # return ac_ref too as varies if N_eng for case 5

            # Push values to output lists
            push!(ac_list, ac)
            Wpointload_fuselage = sum(-point_load.force[3] for point_load in ac.fuselage.point_loads)
            push!(Wpointload_fuselage_list, Wpointload_fuselage)

            # catch e

            # println("An error occurred at iteration $i: ", e)
            # push!(ac_list, NaN)
            # push!(Wpointload_fuselage_list, NaN)
            
            # end

        end

    elseif (case_idx == 31 || case_idx == 631 || case_idx == 32 || case_idx == 632)

        sigma_fc = sigma_fc_min

        # Iterative over independent variable
        for I in CartesianIndices(indep_var_grid)
            i, j = Tuple(I)
            indep_var = indep_var_grid[i, j]

            # try  # FOR DEBUGGING, comment try-catch block

            # Size aircraft
            global ac, ac_ref = PointLoadMethods.analyse_aircraft(
                ac_segment, case_idx, _ac_ref, _Weng_single_ref, indep_var, fcs_loc, N_eng,
                sigma_gt, sigma_fc_wing_min, sigma_fc_wing_max, sigma_fc,  # sigma_fc_min, sigma_fc_max
            )  # return ac_ref too as varies if N_eng for case 5

            # Store values in output grids
            ac_grid[i, j] = ac
            Wpointload_fuselage = sum(-point_load.force[3] for point_load in ac.fuselage.point_loads)
            Wpointload_fuselage_grid[i, j] = Wpointload_fuselage

            # catch e

            # println("An error occurred at iteration $i: ", e)
            # push!(ac_list, NaN)
            # push!(Wpointload_fuselage_list, NaN)
            
            # end

        end

    elseif (case_idx == 4 || case_idx == 64 || case_idx == 5 || case_idx == 65)

        for (i, indep_var_1) in enumerate(indep_var_list_1)

            for (j, indep_var_2) in enumerate(indep_var_list_2)

                sigma_fc = indep_var_1

                # Treat cases separately where one of analyse_aircraft() arguments changes with indep_var
                if (case_idx == 4 || case_idx == 64)
                    fcs_loc = Dict("frac_len" => indep_var_2)
                elseif (case_idx == 5 || case_idx == 65)
                    N_eng = indep_var_2
                end

                # try  # FOR DEBUGGING, comment try-catch block

                # Size aircraft
                global ac, ac_ref = PointLoadMethods.analyse_aircraft(
                    ac_segment, case_idx, _ac_ref, _Weng_single_ref, indep_var_1, fcs_loc, N_eng,
                    sigma_gt, sigma_fc_wing_min, sigma_fc_wing_max, sigma_fc,  #sigma_fc_min, sigma_fc_max
                )  # return ac_ref too as varies if N_eng for case 5

                # Store values in output grids
                ac_grid[i, j] = ac
                Wpointload_fuselage = sum(-point_load.force[3] for point_load in ac.fuselage.point_loads)
                Wpointload_fuselage_grid[i, j] = Wpointload_fuselage

                # catch e

                # println("An error occurred at iteration $i: ", e)
                # push!(ac_list, NaN)
                # push!(Wpointload_fuselage_list, NaN)
                
                # end
            
            end

        end

    end

    result_dict[case_idx]["indep_var_list"] = indep_var_list
    result_dict[case_idx]["indep_var_list_1"] = indep_var_list_1
    result_dict[case_idx]["indep_var_list_2"] = indep_var_list_2
    result_dict[case_idx]["indep_var_grid"] = indep_var_grid
    result_dict[case_idx]["ac_list"] = ac_list
    result_dict[case_idx]["ac_grid"] = ac_grid
    result_dict[case_idx]["ac_ref"] = ac_ref
    result_dict[case_idx]["Wpointload_fuselage"] = Wpointload_fuselage_list

    println("indep_var_list_1: ", indep_var_list_1)
    println("indep_var_list_2: ", indep_var_list_2)

end

# error("message")

## Plotting

using XLSX
using DataFrames
using Tables

xlsx_filename = "exported_results_new.xlsx"

XLSX.openxlsx(xlsx_filename, mode="w") do xf

    # # Remove default empty sheet if present
    # wsnames = XLSX.sheetnames(xf)

    # put this BEFORE the loop over result_dict
    used_default = Ref(false)

    # helper: return existing sheet if present, else:
    # - if Sheet1 exists and hasn't been reused yet, reuse it and rename it via XLSX.rename!
    # - otherwise add a new sheet and return it
    function new_sheet!(xf, sheetname::AbstractString, used_default::Base.RefValue{Bool})
        names = XLSX.sheetnames(xf)

        # if target name already exists, return that worksheet
        if sheetname in names
            idx = findfirst(==(sheetname), names)
            return xf[idx]
        end

        # if default Sheet1 exists and we haven't reused it yet, reuse & rename it
        if !used_default[] && length(names) == 1 && names[1] == "Sheet1"
            ws = xf[1]
            try
                # recommended API to rename a worksheet
                XLSX.rename!(ws, sheetname)
            catch e
                @warn "rename! failed (falling back to reusing Sheet1 without renaming): $e"
                # still reuse the sheet even if rename failed
            end
            used_default[] = true
            return ws
        end

        # otherwise create and return a new sheet
        XLSX.addsheet!(xf, sheetname)
        return xf[XLSX.sheetcount(xf)]
    end

    for (case_idx, case_data) in result_dict
        ac_list = get(case_data, "ac_list", nothing)
        ac_grid = get(case_data, "ac_grid", nothing)
        indep_var_list = get(case_data, "indep_var_list", nothing)
        indep_var_list_1 = collect(get(case_data, "indep_var_list_1", nothing))
        indep_var_list_2 = collect(get(case_data, "indep_var_list_2", nothing))
        indep_var_grid = get(case_data, "indep_var_grid", nothing)
        ac_ref = get(case_data, "ac_ref", nothing)

        # Reference values (skip if missing)
        if ac_ref === nothing || ac_ref === []
            @warn "No ac_ref for case $case_idx, skipping export."
            continue
        end
        WMTO_ref = ac_ref.parg[igWMTO]
        CDS_ref = ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1]

        # # helper: return existing sheet if it exists, else add and return the new sheet
        # function new_sheet!(xf, sheetname::AbstractString)
        #     names = XLSX.sheetnames(xf)
        #     if sheetname in names
        #         idx = findfirst(==(sheetname), names)
        #         return xf[idx]  # reuse existing worksheet
        #     else
        #         XLSX.addsheet!(xf, sheetname)
        #         return xf[XLSX.sheetcount(xf)]
        #     end
        # end

        # 1D case
        if ac_list !== nothing && length(ac_list) > 0 && indep_var_list !== nothing && length(indep_var_list) > 0
            df = DataFrame(
                indep_var = collect(indep_var_list),
                WMTO = [ac.parg[igWMTO] for ac in ac_list],
                CDS = [ac.wing.layout.S * ac.para[iaCD, ipcruise1] for ac in ac_list],
                WMTO_ref = fill(WMTO_ref, length(ac_list)),
                CDS_ref = fill(CDS_ref, length(ac_list)),
            )
            ws = new_sheet!(xf, "case_$(case_idx)_1d", used_default)
            XLSX.writetable!(ws, collect(eachcol(df)), names(df); anchor_cell = XLSX.CellRef("A1"))
        # 2D case
        elseif ac_grid !== nothing && size(ac_grid, 1) > 0 && indep_var_grid !== nothing && size(indep_var_grid, 1) > 0
            n1, n2 = size(ac_grid)
            data = Vector{Dict{Symbol,Any}}()
            for i in 1:n1, j in 1:n2
                ac = ac_grid[i, j]
                indep_var = indep_var_grid[i, j]
                if ac !== nothing && ac !== []
                    push!(data, Dict(
                        :indep_var_1 => (indep_var_list_1 !== nothing && length(indep_var_list_1) >= i) ? indep_var_list_1[i] : "NaN",
                        :indep_var_2 => (indep_var_list_2 !== nothing && length(indep_var_list_2) >= j) ? indep_var_list_2[j] : "NaN",
                        :WMTO => ac.parg[igWMTO],
                        :CDS => ac.wing.layout.S * ac.para[iaCD, ipcruise1],
                        :WMTO_ref => WMTO_ref,
                        :CDS_ref => CDS_ref
                    ))
                end
            end
            df = DataFrame(data)
            ws = new_sheet!(xf, "case_$(case_idx)_2d", used_default)
            XLSX.writetable!(ws, collect(eachcol(df)), names(df); anchor_cell = XLSX.CellRef("A1"))
        end

        # Save independent variable arrays/grids for each case as additional sheets
        if indep_var_list !== nothing && length(indep_var_list) > 0
            ws = new_sheet!(xf, "case_$(case_idx)_indep_var_list", used_default)
            df_ind = DataFrame(indep_var = collect(indep_var_list))
            XLSX.writetable!(ws, collect(eachcol(df_ind)), names(df_ind); anchor_cell = XLSX.CellRef("A1"))
        end

        println("indep_var_list_1: ", indep_var_list_1)
        println("indep_var_list_2: ", indep_var_list_2)
        if indep_var_list_1 !== nothing && length(indep_var_list_1) > 0
            ws = new_sheet!(xf, "case_$(case_idx)_indep_var_list_1", used_default)
            df_ind1 = DataFrame(indep_var_1 = collect(indep_var_list_1))
            XLSX.writetable!(ws, collect(eachcol(df_ind1)), names(df_ind1); anchor_cell = XLSX.CellRef("A1"))
        end

        if indep_var_list_2 !== nothing && length(indep_var_list_2) > 0
            ws = new_sheet!(xf, "case_$(case_idx)_indep_var_list_2", used_default)
            df_ind2 = DataFrame(indep_var_2 = collect(indep_var_list_2))
            XLSX.writetable!(ws, collect(eachcol(df_ind2)), names(df_ind2); anchor_cell = XLSX.CellRef("A1"))
        end

        if indep_var_grid !== nothing && size(indep_var_grid, 1) > 0 && !any(isnan, indep_var_grid)
            println("indep_var_grid: ", indep_var_grid)
            df_grid = DataFrame(indep_var_grid = vec(indep_var_grid))
            ws = new_sheet!(xf, "case_$(case_idx)_indep_var_grid", used_default)
            XLSX.writetable!(ws, collect(eachcol(df_grid)), names(df_grid); anchor_cell = XLSX.CellRef("A1"))
        end
    end
end
