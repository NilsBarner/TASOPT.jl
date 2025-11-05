"""
This script calculates the changes in overall take-off weight, drag, and
ERPK as a result of applying point loads of different magnitude in
different locations and distributions on the airframe. Four variables are
considered:
1) Distrribution of fixed weight (2 kW/kg) along fuselage
2) Distrribution of fixed weight (2 kW/kg) along wing span
3) Variable distrribution of overall-fixed weight (2 kW/kg)
   between wing and fuselage
4) Variable overall weight (specific power 2-4 kW/kg)
"""

using Revise
using TOML
using StaticArrays
using Plots
using TASOPT
include(TASOPT.__TASOPTindices__)

include(joinpath(pwd(), "nils/point_loads/point_load_plots.jl"))
using .PointLoadPlots

include(joinpath(pwd(), "nils/point_loads/system_study/mass/point_load_study_methods_spanfrac.jl"))
using .PointLoadMethods

###

# Top-level settings
ac_segment = "regional"  # "narrowbody" or "regional"
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

_ac_ref, _Weng_single_ref = ac_ref, Weng_single_ref

sigma_fc_range = range(2e3, stop=4e3, length=N_points)
# sigma_fc_range = [3e3]
# N_eng_range = [2, 4, 6, 8]
#####
panel_break_location = 0.37
span_loc_ref = panel_break_location
span_loc_range = collect(range(0.1, stop=1.0, length=N_points))
span_loc_range = sort(unique(vcat(Float64.(span_loc_range), Float64(span_loc_ref))))
# span_loc_range = [span_loc_ref]
#####
fcs_loc_range = range(0.0, stop=1.0, length=N_points)
# fcs_loc_range = [1.0]
wing_frac_range = range(0.0, stop=1.0, length=N_points)
# wing_frac_range = [1.0]

indep_var_grid = Array{Any}(
    # undef, length(sigma_fc_range), length(N_eng_range), length(fcs_loc_range), length(wing_frac_range)
    undef, length(sigma_fc_range), length(span_loc_range), length(fcs_loc_range), length(wing_frac_range)
)
ac_grid = Array{Any}(
    # undef, length(sigma_fc_range), length(N_eng_range), length(fcs_loc_range), length(wing_frac_range)
    undef, length(sigma_fc_range), length(span_loc_range), length(fcs_loc_range), length(wing_frac_range)
)
Wpointload_fuselage_grid = fill(
    # NaN, length(sigma_fc_range), length(N_eng_range), length(fcs_loc_range), length(wing_frac_range)
    NaN, length(sigma_fc_range), length(span_loc_range), length(fcs_loc_range), length(wing_frac_range)
)

for (i, sigma_fc) in enumerate(sigma_fc_range)
    # for (j, N_eng) in enumerate(N_eng_range)
    for (j, span_loc) in enumerate(span_loc_range)
        for (k, fcs_loc) in enumerate(fcs_loc_range)
            for (l, wing_frac) in enumerate(wing_frac_range)
                # indep_var_grid[i, j, k, l] = (sigma_fc, N_eng, fcs_loc, wing_frac)
                indep_var_grid[i, j, k, l] = (sigma_fc, span_loc, fcs_loc, wing_frac)

                _fcs_loc = Dict("frac_len" => fcs_loc)
                _span_loc = Dict("frac_span" => span_loc)

                try  # FOR DEBUGGING, comment try-catch block

                    # Size aircraft
                    global ac, ac_ref = PointLoadMethods.analyse_aircraft(
                        # ac_segment, _ac_ref, _Weng_single_ref, _fcs_loc, N_eng,
                        ac_segment, _ac_ref, _Weng_single_ref, _fcs_loc, _span_loc,
                        sigma_gt, sigma_fc_wing_max, sigma_fc, wing_frac,
                    )

                    # Store values in output grids
                    ac_grid[i, j, k, l] = ac
                    Wpointload_fuselage = sum(-point_load.force[3] for point_load in ac.fuselage.point_loads)
                    Wpointload_fuselage_grid[i, j, k, l] = Wpointload_fuselage
                    # println("Wpointload_fuselage: ", Wpointload_fuselage)

                catch e

                    # println("sigma_fc: $sigma_fc, N_eng: $N_eng, fcs_loc: $fcs_loc, wing_frac: $wing_frac")
                    println("sigma_fc: $sigma_fc, _span_loc: $_span_loc, fcs_loc: $fcs_loc, wing_frac: $wing_frac")

                end

            end

        end
    
    end

end

error("message")  # NILS: stop here for debugging

### Save study results to .xlsx file for plotting in python

using XLSX
using DataFrames
using Tables

xlsx_filename = "point_load_study_results_regional_spanfrac.xlsx"

# Prepare vectors to collect flattened data
indep_var_flat = Vector{Any}()
CDS_flat = Vector{Any}()
WMTO_flat = Vector{Any}()
CDS_ref_flat = Vector{Any}()
WMTO_ref_flat = Vector{Any}()
Wpointload_fuselage_flat = Vector{Any}()

n1, n2, n3, n4 = size(indep_var_grid)
for i in 1:n1, j in 1:n2, k in 1:n3, l in 1:n4
    indep = indep_var_grid[i, j, k, l]
    if isassigned(ac_grid, i, j, k, l)
        ac = ac_grid[i, j, k, l]
        if ac !== nothing && ac !== []
            push!(indep_var_flat, indep)
            push!(CDS_flat, ac.wing.layout.S * ac.para[iaCD, ipcruise1])
            push!(WMTO_flat, ac.parg[igWMTO])
            push!(CDS_ref_flat, ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1])
            push!(WMTO_ref_flat, ac_ref.parg[igWMTO])
            push!(Wpointload_fuselage_flat, Wpointload_fuselage_grid[i, j, k, l])
            continue
        end
    end
    # If not assigned or invalid, push NaN
    push!(indep_var_flat, indep)
    push!(CDS_flat, "NaN")
    push!(WMTO_flat, "NaN")
    push!(CDS_ref_flat, "NaN")
    push!(WMTO_ref_flat, "NaN")
    push!(Wpointload_fuselage_flat, "NaN")
end

# Split indep_var_flat into columns
sigma_fc_col = [x[1] for x in indep_var_flat]
# N_eng_col = [x[2] for x in indep_var_flat]
span_loc_col = [x[2] for x in indep_var_flat]
fcs_loc_col = [x[3] for x in indep_var_flat]
wing_frac_col = [x[4] for x in indep_var_flat]

df = DataFrame(
    sigma_fc = sigma_fc_col,
    # N_eng = N_eng_col,
    span_loc = span_loc_col,
    fcs_loc = fcs_loc_col,
    wing_frac = wing_frac_col,
    CDS = CDS_flat,
    CDS_ref = CDS_ref_flat,
    WMTO = WMTO_flat,
    WMTO_ref = WMTO_ref_flat,
    Wpointload_fuselage = Wpointload_fuselage_flat
)

# Write DataFrame to Excel using the filename-based API (correct signature)
XLSX.writetable(
    xlsx_filename,
    Tables.columntable(df);
    sheetname = "results",
    overwrite = true
)

### Plot aircraft top-views for different cases (uncomment and run one-by-one)

using Colors

ac_ref = ac_grid[1, 1, 1, length(wing_frac_range)]
ac_1 = ac_grid[1, 1, :, 1]
ac_2 = ac_grid[1, length(N_eng_range), 1, length(wing_frac_range)]
ac_3 = ac_grid[1, 1, length(fcs_loc_range), :]
ac_4 = ac_grid[length(sigma_fc_range), 1, 1, length(wing_frac_range)]

include(joinpath(pwd(), "nils/point_loads/point_load_plots.jl"))  # Adjust the path as needed
using .PointLoadPlots  # Import the module

# ac_list_flat, colour_list = [ac_ref, ac_1[1], ac_1[4]], [colorant"#206095", colorant"#206095"]
# ac_list_flat, colour_list = [ac_ref, ac_2], [colorant"#a8bd3a"]  # if want to plot all cases
# ac_list_flat, colour_list = [ac_ref, ac_3[1], ac_3[2], ac_3[3], ac_3[4]], [colorant"#871a5b", colorant"#871a5b", colorant"#871a5b", colorant"#871a5b"]  # if want to plot all cases
ac_list_flat, colour_list = [ac_ref, ac_4], [colorant"#f66068"]  # if want to plot all cases

p = PointLoadPlots.compare_stickfig(
    Vector{Any}(ac_list_flat); colour_list=colour_list, annotate_text = false, annotate_group = false, framestyle = :none,
)
# savefig(p, "stickfig_case_1.svg")
# savefig(p, "stickfig_case_2.svg")
# savefig(p, "stickfig_case_3.svg")
savefig(p, "stickfig_case_4.svg")
