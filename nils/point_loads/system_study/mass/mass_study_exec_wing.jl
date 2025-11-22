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
using DataFrames
using CSV
using TASOPT
include(TASOPT.__TASOPTindices__)

include(joinpath(pwd(), "nils", "point_loads", "point_load_plots.jl"))
using .PointLoadPlots

include(joinpath(pwd(), "nils", "point_loads", "system_study", "mass", "mass_study_methods_wing.jl"))
using .PointLoadMethods

###

# Top-level settings
ac_segment = "narrowbody"  # "narrowbody" or "regional"
N_points = 20

###

# Model reference aircraft
ac_ref, Weng_single_ref = PointLoadMethods.analyse_reference_aircraft(ac_segment)
# ac_ref, Weng_single_ref = nothing, nothing

# Calculate maximum propulsive power throughout mission
T_net_max_ref = ac_ref.pare[ieFe, :] * ac_ref.parg[igneng]  # NILS: see calculate_thrust_from_ROC!() for proof that pare[ieFe] stands for per-engine thrust
V_0_max_ref = ac_ref.pare[ieu0, :]
P_prop_max_ref = maximum(T_net_max_ref .* V_0_max_ref)

df = DataFrame(
    CDS = ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1],
    WMTO = ac_ref.parg[igWMTO],
    Vol_wing = 2*(ac_ref.wing.center.volume + ac_ref.wing.inboard.volume + ac_ref.wing.outboard.volume),
    PFEI = ac_ref.parm[imPFEI],
    seats_abreast = ac_ref.fuselage.cabin.seats_abreast_main,
    Vol_nacelle = ac_ref.nils.V_fcs_nacelle,
    L_fuse = ac_ref.fuselage.layout.l_cabin_cylinder,
    span = ac_ref.wing.layout.span,
    Vol_fuse = ac_ref.nils.V_fcs_fuselage,
    length = ac_ref.fuselage.layout.x_end,
    P_prop_max = P_prop_max_ref,

    # root_span = root_span_vec,
    # dfan = dfan_vec,

    Wfuse = ac_ref.fuselage.weight,
    Wwing = ac_ref.wing.weight,
    Whtail = ac_ref.htail.weight,
    Wvtail = ac_ref.vtail.weight,
    Wfuseshell = ac_ref.fuselage.shell.weight.W,
    Wfusebendv = ac_ref.fuselage.bendingmaterial_v.weight.W,
    Wfusebendh = ac_ref.fuselage.bendingmaterial_h.weight.W,
    Wwingcentre = ac_ref.wing.center.weight,
    Wwinginboard = ac_ref.wing.inboard.weight,
    Wwingoutboard = ac_ref.wing.outboard.weight,

    CLS = ac_ref.wing.layout.S * ac_ref.para[iaCL, ipcruise1],
    CDfuseS = ac_ref.wing.layout.S * ac_ref.para[iaCDfuse, ipcruise1],
    CDiS = ac_ref.wing.layout.S * ac_ref.para[iaCDi, ipcruise1],
    CDwingS = ac_ref.wing.layout.S * ac_ref.para[iaCDwing, ipcruise1],
    CDhtailS = ac_ref.wing.layout.S * ac_ref.para[iaCDhtail, ipcruise1],
    CDvtailS = ac_ref.wing.layout.S * ac_ref.para[iaCDvtail, ipcruise1],

    # eta_grav_tank = etagravtank_vec,
)

CSV.write(joinpath(pwd(), "nils", "point_loads", "system_study", "mass", "data", "point_load_study_results_211125_ref.csv"), df)
# error("stop here for debugging")  # NILS

###

function run_study(study_idx::Int)

    # result vectors
    study_idx_vec = Int[]
    index_vec = Vector{NTuple{5,Int}}()
    sigma_fcs_vec = Float64[]
    span_loc_vec = Any[]
    fcs_loc_vec = Any[]
    wing_frac_vec = Float64[]
    nacelle_frac_vec = Float64[]

    CDS_vec = Float64[]
    WMTO_vec = Float64[]
    Vol_wing_vec = Float64[]
    PFEI_vec = Float64[]
    seats_abreast_vec = Int[]
    Vol_nacelle_vec = Float64[]
    L_fuse_vec = Float64[]
    span_vec = Float64[]
    Vol_fuse_vec = Float64[]
    length_vec = Float64[]
    P_prop_max_vec = Float64[]

    Wfuse_vec = Float64[]
    Wwing_vec = Float64[]
    Whtail_vec = Float64[]
    Wvtail_vec = Float64[]
    Wfuseshell_vec = Float64[]
    Wfusebendv_vec = Float64[]
    Wfusebendh_vec = Float64[]
    Wwingcentre_vec = Float64[]
    Wwinginboard_vec = Float64[]
    Wwingoutboard_vec = Float64[]

    CLS_vec = Float64[]
    CDfuseS_vec = Float64[]
    CDiS_vec = Float64[]
    CDwingS_vec = Float64[]
    CDhtailS_vec = Float64[]
    CDvtailS_vec = Float64[]

    if study_idx == 1
        # sigma_fcs_range = [3.25e3]
        sigma_fcs_range = [4e3]
        span_loc_range = [-1.0]
        fcs_loc_range = ["eng"]
        wing_frac_range = [0.0]
        nacelle_frac_range = range(0.2, stop=1.0, length=N_points)
    elseif study_idx == 2
        # sigma_fcs_range = [3.25e3]
        sigma_fcs_range = [4e3]
        span_loc_range = ["eng"]
        fcs_loc_range = [-1.0]
        wing_frac_range = [1.0]
        nacelle_frac_range = range(0.2, stop=1.0, length=N_points)
    elseif study_idx == 3
        # sigma_fcs_range = [3.25e3]
        sigma_fcs_range = [4e3]
        span_loc_range = [-1.0]
        fcs_loc_range = range(0.01, stop=0.99, length=N_points)
        wing_frac_range = [0.0]
        nacelle_frac_range = [0.2]
    elseif study_idx == 4
        # sigma_fcs_range = [3.25e3]
        sigma_fcs_range = [4e3]
        span_loc_range = range(0.01, stop=0.99, length=N_points)
        fcs_loc_range = [-1.0]
        wing_frac_range = [1.0]
        nacelle_frac_range = [0.2]
    elseif study_idx == 5
        sigma_fcs_range = [4e3]
        span_loc_range = range(0.01, stop=0.99, length=N_points)
        fcs_loc_range = range(0.01, stop=0.99, length=N_points)
        # wing_frac_range = range(0.01, stop=0.99, length=N_points)
        wing_frac_range = [0.5]
        nacelle_frac_range = [0.2]
    elseif study_idx == 6
        sigma_fcs_range = range(1.5e3, stop=4e3, length=N_points)
        span_loc_range = range(0.01, stop=0.99, length=N_points)
        fcs_loc_range = [-1.0]
        wing_frac_range = [1.0]
        nacelle_frac_range = [0.2]
    elseif study_idx == 7
        sigma_fcs_range = range(1.5e3, stop=4e3, length=N_points)
        span_loc_range = [-1.0]
        fcs_loc_range = range(0.01, stop=0.99, length=N_points)
        wing_frac_range = [0.0]
        nacelle_frac_range = [0.2]
    elseif study_idx == 8
        sigma_fcs_range = range(1.5e3, stop=4e3, length=N_points)
        span_loc_range = [-1.0]
        fcs_loc_range = [-1.0]
        wing_frac_range = [-1.0]
        nacelle_frac_range = [1.0]
    end

    ###

    indep_var_grid = Array{Any}(
        undef, length(sigma_fcs_range), length(span_loc_range), length(fcs_loc_range), length(wing_frac_range), length(nacelle_frac_range),
    )
    ac_grid = Array{Any}(
        undef, length(sigma_fcs_range), length(span_loc_range), length(fcs_loc_range), length(wing_frac_range), length(nacelle_frac_range),
    )
    Wpointload_fuselage_grid = fill(
        NaN, length(sigma_fcs_range), length(span_loc_range), length(fcs_loc_range), length(wing_frac_range), length(nacelle_frac_range),
    )
    Wpointload_wing_grid = fill(
        NaN, length(sigma_fcs_range), length(span_loc_range), length(fcs_loc_range), length(wing_frac_range), length(nacelle_frac_range),
    )
    Wpointload_nacelle_grid = fill(
        NaN, length(sigma_fcs_range), length(span_loc_range), length(fcs_loc_range), length(wing_frac_range), length(nacelle_frac_range),
    )

    for (i, sigma_fcs) in enumerate(sigma_fcs_range)
        for (j, span_loc) in enumerate(span_loc_range)
            for (k, fcs_loc) in enumerate(fcs_loc_range)
                for (l, wing_frac) in enumerate(wing_frac_range)
                    for (m, nacelle_frac) in enumerate(nacelle_frac_range)
                        indep_var_grid[i, j, k, l, m] = (sigma_fcs, span_loc, fcs_loc, wing_frac, nacelle_frac)

                        try  # FOR DEBUGGING, comment try-catch block

                            # Size aircraft
                            global ac, ac_ref = PointLoadMethods.analyse_aircraft(
                                ac_segment, ac_ref, Weng_single_ref,
                                fcs_loc, span_loc, sigma_fcs, wing_frac, nacelle_frac,
                            )

                            # # Store values in output grids
                            # ac_grid[i, j, k, l] = ac
                            # Wpointload_fuselage = sum(-point_load.force[3] for point_load in ac.fuselage.point_loads)  # total in fuselage
                            # Wpointload_fuselage_grid[i, j, k, l] = Wpointload_fuselage
                            # Wpointload_wing = sum(-point_load.force[3] for point_load in ac.wing.point_loads) * 2  # both wing halfs
                            # Wpointload_wing_grid[i, j, k, l] = Wpointload_wing
                            # # Wpointload_nacelle = -ac.engine.model.point_load * ac.parg[igneng]  # all nacelles
                            # Wpointload_nacelle = -ac.engine.point_load * ac.parg[igneng]  # all nacelles
                            # Wpointload_nacelle_grid[i, j, k, l] = Wpointload_nacelle

                            # Calculate maximum propulsive power throughout mission
                            T_net_max = ac.pare[ieFe, :] * ac.parg[igneng]  # NILS: see calculate_thrust_from_ROC!() for proof that pare[ieFe] stands for per-engine thrust
                            V_0_max = ac.pare[ieu0, :]
                            P_prop_max = maximum(T_net_max .* V_0_max)

                            # save results immediately
                            push!(study_idx_vec, study_idx)
                            push!(index_vec, (i,j,k,l,m))
                            push!(sigma_fcs_vec, sigma_fcs)
                            push!(span_loc_vec, span_loc)
                            push!(fcs_loc_vec, fcs_loc)
                            push!(wing_frac_vec, wing_frac)
                            push!(nacelle_frac_vec, nacelle_frac)

                            push!(CDS_vec, ac.wing.layout.S * ac.para[iaCD, ipcruise1])
                            push!(WMTO_vec, ac.parg[igWMTO])
                            push!(Vol_wing_vec, 2*(ac.wing.center.volume + ac.wing.inboard.volume + ac.wing.outboard.volume))
                            push!(PFEI_vec, ac.parm[imPFEI])
                            push!(seats_abreast_vec, ac.fuselage.cabin.seats_abreast_main)
                            push!(Vol_nacelle_vec, ac.nils.V_fcs_nacelle)  # ac.parg[igVfcsavnacetot])
                            push!(L_fuse_vec, ac.fuselage.layout.l_cabin_cylinder)
                            push!(span_vec, ac.wing.layout.span)
                            push!(Vol_fuse_vec, ac.nils.V_fcs_fuselage)  # ac.parg[igVfcsfus])
                            push!(length_vec, ac.fuselage.layout.x_end)
                            push!(P_prop_max_vec, P_prop_max)

                            push!(Wfuse_vec, ac.fuselage.weight)
                            push!(Wwing_vec, ac.wing.weight)
                            push!(Whtail_vec, ac.htail.weight)
                            push!(Wvtail_vec, ac.vtail.weight)
                            push!(Wfuseshell_vec, ac.fuselage.shell.weight.W)
                            push!(Wfusebendv_vec, ac.fuselage.bendingmaterial_v.weight.W)
                            push!(Wfusebendh_vec, ac.fuselage.bendingmaterial_h.weight.W)
                            push!(Wwingcentre_vec, ac.wing.center.weight)
                            push!(Wwinginboard_vec, ac.wing.inboard.weight)
                            push!(Wwingoutboard_vec, ac.wing.outboard.weight)

                            push!(CLS_vec, ac.wing.layout.S * ac.para[iaCL, ipcruise1])
                            push!(CDfuseS_vec, ac.wing.layout.S * ac.para[iaCDfuse, ipcruise1])
                            push!(CDiS_vec, ac.wing.layout.S * ac.para[iaCDi, ipcruise1])
                            push!(CDwingS_vec, ac.wing.layout.S * ac.para[iaCDwing, ipcruise1])
                            push!(CDhtailS_vec, ac.wing.layout.S * ac.para[iaCDhtail, ipcruise1])
                            push!(CDvtailS_vec, ac.wing.layout.S * ac.para[iaCDvtail, ipcruise1])

                        catch e

                            println(e)
                            # println("sigma_fcs: $sigma_fcs, _span_loc: $_span_loc, fcs_loc: $fcs_loc, wing_frac: $wing_frac")

                        end

                        # error("stop here for debugging")  # NILS: stop here for debugging

                    end
                end
            end
        end
    end

    df = DataFrame(
        study_idx = study_idx_vec,
        index = index_vec,
        sigma_fcs = sigma_fcs_vec,
        span_loc = span_loc_vec,
        fcs_loc = fcs_loc_vec,
        wing_frac = wing_frac_vec,
        nacelle_frac = nacelle_frac_vec,

        CDS = CDS_vec,
        WMTO = WMTO_vec,
        Vol_wing = Vol_wing_vec,
        PFEI = PFEI_vec,
        seats_abreast = seats_abreast_vec,
        Vol_nacelle = Vol_nacelle_vec,
        L_fuse = L_fuse_vec,
        span = span_vec,
        Vol_fuse = Vol_fuse_vec,
        length = length_vec,
        P_prop_max = P_prop_max_vec,

        # root_span = root_span_vec,
        # dfan = dfan_vec,

        Wfuse = Wfuse_vec,
        Wwing = Wwing_vec,
        Whtail = Whtail_vec,
        Wvtail = Wvtail_vec,
        Wfuseshell = Wfuseshell_vec,
        Wfusebendv = Wfusebendv_vec,
        Wfusebendh = Wfusebendh_vec,
        Wwingcentre = Wwingcentre_vec,
        Wwinginboard = Wwinginboard_vec,
        Wwingoutboard = Wwingoutboard_vec,

        CLS = CLS_vec,
        CDfuseS = CDfuseS_vec,
        CDiS = CDiS_vec,
        CDwingS = CDwingS_vec,
        CDhtailS = CDhtailS_vec,
        CDvtailS = CDvtailS_vec,

        # eta_grav_tank = etagravtank_vec,
    )

    CSV.write(joinpath(pwd(), "nils", "point_loads", "system_study", "mass", "data", "point_load_study_results_211125_$(study_idx).csv"), df)

end

for study_idx in 1:8
    run_study(study_idx)
end

# using Plots
# p = TASOPT.stickfig(ac)  # line commented by Nils
# display(p)
error("message")  # NILS: stop here for debugging



### Save study results to .xlsx file for plotting in python

using XLSX
using DataFrames
using Tables

xlsx_filename = "point_load_study_results_regional_spanfrac_wing_new.xlsx"

# Prepare vectors to collect flattened data
indep_var_flat = Vector{Any}()
CDS_flat = Vector{Any}()
WMTO_flat = Vector{Any}()
CDS_ref_flat = Vector{Any}()
WMTO_ref_flat = Vector{Any}()
Wpointload_fuselage_flat = Vector{Any}()
Wpointload_wing_flat = Vector{Any}()
Wpointload_nacelle_flat = Vector{Any}()
Wwing_flat = Vector{Any}()

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
            push!(Wpointload_wing_flat, Wpointload_wing_grid[i, j, k, l])
            push!(Wpointload_nacelle_flat, Wpointload_nacelle_grid[i, j, k, l])
            push!(Wwing_flat, ac.wing.weight)
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
    push!(Wpointload_wing_flat, "NaN")
    push!(Wpointload_nacelle_flat, "NaN")
    push!(Wwing_flat, "NaN")
end

# Split indep_var_flat into columns
sigma_fcs_col = [x[1] for x in indep_var_flat]
span_loc_col = [x[2] for x in indep_var_flat]
fcs_loc_col = [x[3] for x in indep_var_flat]
wing_frac_col = [x[4] for x in indep_var_flat]

df = DataFrame(
    sigma_fcs = sigma_fcs_col,
    # N_eng = N_eng_col,
    span_loc = span_loc_col,
    fcs_loc = fcs_loc_col,
    wing_frac = wing_frac_col,
    CDS = CDS_flat,
    CDS_ref = CDS_ref_flat,
    WMTO = WMTO_flat,
    WMTO_ref = WMTO_ref_flat,
    Wpointload_fuselage = Wpointload_fuselage_flat,
    Wpointload_wing = Wpointload_wing_flat,
    Wpointload_nacelle = Wpointload_nacelle_flat,
    Wwing = Wwing_flat,
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
ac_4 = ac_grid[length(sigma_fcs_range), 1, 1, length(wing_frac_range)]

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
