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

# Top-level settings
ac_segment = "narrowbody"  # "narrowbody" or "regional"
N_points = 2

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
case_idx_list = [5]
# case_idx_list = [1, 61]
# case_idx_list = [1, 21, 22, 31, 32, 4, 5, 61, 621, 622, 631, 632, 64, 65]

result_dict = Dict(
    i => Dict(
        "ac" => [],
        "ac_ref" => Union{Nothing, aircraft},
        "Wpointload_fuselage" => [],
    ) for i in case_idx_list
)

###

# Iterate over cases
for (h, case_idx) in enumerate(case_idx_list)

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
    elseif (case_idx == 21 || case_idx == 621)
        indep_var_list = range(W_fc_min, stop=W_fc_max, length=N_points)  # entire FCS weight varies
        fcs_loc = "apu"
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
    elseif (case_idx == 22 || case_idx == 622)
        indep_var_list = range(W_fc_min, stop=W_fc_max, length=N_points)  # entire FCS weight varies
        fcs_loc = "eng"
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
    elseif (case_idx == 31 || case_idx == 631)
        indep_var_list = range(W_fc_wing_min, stop=W_fc_wing_max, length=N_points)  # fraction of FCS weight in nacelle varies, and thus inversely that in fuselage for fixed overall weight
        fcs_loc = "apu"
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
    elseif (case_idx == 32 || case_idx == 632)
        indep_var_list = range(W_fc_wing_min, stop=W_fc_wing_max, length=N_points)  # fraction of FCS weight in nacelle varies, and thus inversely that in fuselage for fixed overall weight
        fcs_loc = "eng"
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
    elseif (case_idx == 4 || case_idx == 64)
        indep_var_list = range(0, stop=1, length=N_points)  # location of point load varies
        fcs_loc = 0.0  # value not of importance
        _ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref
        N_eng = 2
    elseif (case_idx == 5 || case_idx == 65)
        indep_var_list = [2, 4, 6, 8]
        fcs_loc = 0.0  # value not of importance
        _ac_ref, _Weng_single_ref = nothing, nothing
        N_eng = 0  # value not of importance
    end

    # Output lists
    global ac_list = []
    global Wpointload_fuselage_list = Float64[]

    # Iterative over independent variable
    for (i, indep_var) in enumerate(indep_var_list)

        # Treat cases separately where one of analyse_aircraft() arguments changes with indep_var
        if (case_idx == 4 || case_idx == 64)
            fcs_loc = Dict("frac_len" => indep_var)
        elseif (case_idx == 5 || case_idx == 65)
            N_eng = indep_var
        end

        # try  # FOR DEBUGGING, comment try-catch block

        # Size aircraft
        sigma_fc = sigma_fc_min  # consider heaviest FCS when sigma_fc not indep_var
        global ac, ac_ref = PointLoadMethods.analyse_aircraft(
            ac_segment, case_idx, _ac_ref, _Weng_single_ref, indep_var, fcs_loc, N_eng,
            sigma_gt, sigma_fc_wing_min, sigma_fc_wing_max, sigma_fc,
        )  # return ac_ref too as varies if N_eng for case 5

        # Push values to output lists
        push!(ac_list, ac)
        Wpointload_fuselage = sum(-point_load.force[3] for point_load in ac.fuselage.point_loads)
        push!(Wpointload_fuselage_list, Wpointload_fuselage)
        # println("Wpointload_fuselage: ", Wpointload_fuselage)  # NILS: this was added on 20.10.2025 to check correct implementation of fuselage point load

        # catch e

        # println("An error occurred at iteration $i: ", e)
        # push!(ac_list, NaN)
        # push!(Wpointload_fuselage_list, NaN)
        
        # end

    end

    result_dict[case_idx]["ac"] = ac_list
    result_dict[case_idx]["ac_ref"] = ac_ref
    result_dict[case_idx]["Wpointload_fuselage"] = Wpointload_fuselage_list

end

# error("message")  # NILS: stop here for debugging

### Plotting

# Stickfigure plots

include(joinpath(pwd(), "nils/point_loads/point_load_plots.jl"))  # Adjust the path as needed
using .PointLoadPlots  # Import the module

# colour_list = [:grey, :grey, :blue, :blue, :green, :green, :purple, :purple, :red, :red, :orange, :orange, :brown, :brown]
colour_list = [:red]  # if want to plot only two aircraft at a time (one case)
ac_list_flat = vcat([[result_dict[i]["ac"][1], result_dict[i]["ac"][end]] for i in case_idx_list]...)  # if want to plot all cases

# p = PointLoadPlots.compare_stickfig(
#     Vector{Any}(ac_list_flat); colour_list=colour_list, annotate_text = false, annotate_group = false
# )
# savefig(p, "stickfig_case_$(string(case_idx_list[1])).png")
# display(p)

# error("message")

###

# Weight for different segments

# ac_refs_flat = vcat([result_dict[i]["ac_ref"] for i in case_idx_list]...)

# PointLoadPlots.plot_wrel_vs_segment(
#     Vector{Any}(ac_refs_flat), result_dict, [1, 61], N_points,
#     sigma_fc_min, sigma_fc_max,
# )

# error("message")

###

# Weight for different fuel cell system specific powers

# PointLoadPlots.plot_wrel_vs_sigma_fcs(
#     ac_ref, result_dict, [1, 21, 22], N_points,
#     sigma_fc_min, sigma_fc_max,
# )

# PointLoadPlots.plot_wrel_vs_sigma_fcs(
#     ac_ref, result_dict, [31, 32], N_points,
#     sigma_fc_min, sigma_fc_max,
# )

# PointLoadPlots.plot_wmto_vs_sigma_fcs(
#     ac_ref, result_dict, N_points,
#     sigma_fc_min, sigma_fc_max,
# )

# error("message")

# ###

# Weight for different weight splits

# PointLoadPlots.plot_weight_drag_vs_weight_split(
#     ac_ref, result_dict, N_points,
#     Weng_single_ref, sigma_gt, sigma_fc_wing_min, sigma_fc_wing_max,
# )

# error("message")

# ###

# Weight for different weight positions (fuselage and wing)

# PointLoadPlots.plot_wmto_vs_weight_position(
#     ac_ref, result_dict, N_points, length([2, 4, 6, 8]),
# )

# error("message")

###

# Drag for different fuel cell system specific powers

# ac_refs_flat = vcat([result_dict[i]["ac_ref"] for i in case_idx_list]...)

# PointLoadPlots.plot_drag_vs_sigma_fcs(
#     Vector{Any}(ac_refs_flat), result_dict, [1, 61], N_points,
#     sigma_fc_min, sigma_fc_max,
# )

###

# Weight and drag for different segments

ac_refs_flat = vcat([result_dict[i]["ac_ref"] for i in case_idx_list]...)

PointLoadPlots.scatter_weight_drag_segment(
    Vector{Any}(ac_refs_flat), result_dict, case_idx_list,
)