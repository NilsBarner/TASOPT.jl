using Revise
using TOML
using StaticArrays
using Plots
using TASOPT
include(TASOPT.__TASOPTindices__)

include(joinpath(pwd(), "nils/point_loads/point_load_plots.jl"))
using .PointLoadPlots

# include(joinpath(pwd(), "nils/point_loads/point_load_study_methods.jl"))
# using .PointLoadMethods

###

# Top-level settings
ac_segment = "narrowbody"  # "narrowbody" or "regional"
##### NILS: RESET "REGIONAL" - CURRENTLY DOES NOT RUN, PRESUMABLY BECAUSE OF INFEASIBLE INPUT FROM MASS-EXPERIMENTS #####
N_points = 2
fcs_loc = "nacelle"  # "nacelle", "wing", "fuselage"

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
    # toml_data["Propulsion"]["V_fcs"] = V_fcs
    toml_data["Propulsion"]["V_fcs_nacelle"] = 0.0
    toml_data["Wing"]["V_fcs_wing"] = 0.0
    toml_data["Fuselage"]["V_fcs_fuselage"] = 0.0
    toml_data["Fuselage"]["fcs_fuselage_location"] = 1.0  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl"
    toml_data["Wing"]["AR"] = 10.1
    toml_data["Wing"]["tdivc_scale"] = 1.0
    toml_data["Wing"]["TR_scale"] = 1.0
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

    return ac_ref, Weng_single_ref

end

###

# Model reference aircraft
ac_ref, Weng_single_ref = analyse_reference_aircraft(ac_segment)

###

P_fcs = 20e6  # [W] total FCS power

fcs_locs = ["nacelle", "wing", "fuselage"]
radius_range = range(90, 140, length=N_points) * 0.0254  # [m]
AR_range = range(6, 18, length=N_points)
TR_scale_range = range(0.5, 2, length=N_points)
tdivc_scale_range = range(0.75, 1.5, length=N_points)
Neng_range = [2, 4, 6, 8]
Vspec_range = range(3, 5.5, length=N_points) * 1e3 * 1e3  # [m^3]  Table 2 in FZO-PPN-REP-0031 - Fuel Cells Technical Report (STACK)
# Vspec_range = range(0.2, 0.4, length=N_points) * 1e3 * 1e3  # [m^3]  Table 2 in FZO-PPN-REP-0031 - Fuel Cells Technical Report (SYSTEM)
Vspec_range = range(0.2, 1.0, length=N_points) * 1e3 * 1e3  # NILS' estimate
HTRf_range = range(0.2, 0.9, length=N_points)

# indep_var_grid = Array{Any}(
#     undef, length(fcs_locs), length(radius_range), length(AR_range), length(Neng_range), length(Vspec_range)
# )
# ac_grid = Array{Any}(
#     undef, length(fcs_locs), length(radius_range), length(AR_range), length(Neng_range), length(Vspec_range)
# )

# for (h, fcs_loc) in enumerate(fcs_locs)
#     for (i, radius) in enumerate(radius_range)
#         for (j, AR) in enumerate(AR_range)
#             for (k, N_eng) in enumerate(Neng_range)
#                 for (l, Vspec) in enumerate(Vspec_range)

#                     indep_var_grid[h, i, j, k, l] = (fcs_loc, radius, AR, N_eng, Vspec)

fcs_loc = "nacelle"
fcs_fuselage_location = 0.0
radius = 100 * 0.0254  # [m]
AR = 10.1
N_eng = 2
Vspec = 100e6
TR_scale = 1.0
tdivc_scale = 1.0

ac_list = []
Vol_wing_list = []
Vol_nace_list = []

# Vspec_range = [0.5e6]

# for Vspec in Vspec_range
# for AR in AR_range
# for TR_scale in TR_scale_range
# for tdivc_scale in tdivc_scale_range
# for N_eng in Neng_range
for HTR_f in HTRf_range

    global V_fcs = P_fcs / Vspec  # [m^3] total FCS volume

    ###

    # Write engine point loads to .toml file
    if ac_segment == "narrowbody"
        toml_data = TOML.parsefile(joinpath(pwd(), "example", "cryo_input.toml"))
    elseif ac_segment == "regional"
        toml_data = TOML.parsefile(joinpath(pwd(), "nils", "point_loads", "default_regional_cryo.toml"))
    end
    toml_data["Propulsion"]["number_of_engines"] = N_eng
    toml_data["Propulsion"]["Turbomachinery"]["HTR_fan"] = HTR_f
    # toml_data["Propulsion"]["V_fcs"] = V_fcs
    if fcs_loc == "nacelle"
        toml_data["Propulsion"]["V_fcs_nacelle"] = V_fcs / N_eng
        toml_data["Propulsion"]["V_fcs_wing"] = 0.0
        toml_data["Propulsion"]["V_fcs_fuselage"] = 0.0
    elseif fcs_loc == "wing"
        toml_data["Wing"]["V_fcs_nacelle"] = 0.0
        toml_data["Wing"]["V_fcs_wing"] = V_fcs
        toml_data["Wing"]["V_fcs_fuselage"] = 0.0
    elseif fcs_loc == "fuselage"
        toml_data["Fuselage"]["V_fcs_nacelle"] = 0.0
        toml_data["Fuselage"]["V_fcs_wing"] = 0.0
        toml_data["Fuselage"]["V_fcs_fuselage"] = V_fcs
        toml_data["Fuselage"]["fcs_fuselage_location"] = fcs_fuselage_location  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl"
    end
    toml_data["Wing"]["AR"] = AR
    toml_data["Wing"]["tdivc_scale"] = tdivc_scale
    toml_data["Wing"]["TR_scale"] = TR_scale
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

    # ###############################################################################################
    # # Add fuselage point loads to mutable structure
    # _fcs_loc = Dict("frac_span" => 0.5)
    # Fz_point_wing = 0.0
    # wing_load = TASOPT.structures.PointLoad(
    #     force = SVector(0.0, 0.0, Fz_point_wing),  # NILS: previously, I mistakenly multiplied by number of engines here (20.10.2025)
    #     r = SVector(0.0, _fcs_loc, 0.0),
    #     frame = TASOPT.WORLD
    # )
    # TASOPT.structures.add_wing_point_load!(ac.wing, wing_load)
    # ###############################################################################################

    try  # FOR DEBUGGING, comment try-catch block

        # Size aircraft
        size_aircraft!(ac, iter=100, printiter=false)  # disabled default iteration print-outs to console
        # # Store aircraft object in output grids
        # ac_grid[h, i, j, k, l] = ac

        # Push values to output lists
        push!(ac_list, ac)
        push!(Vol_wing_list, 2 * (ac.wing.center.volume + ac.wing.inboard.volume + ac.wing.outboard.volume))
        push!(Vol_nace_list, ac.parg[igVfcsavnacetot])

    catch e

        println("fcs_loc: $fcs_loc, radius: $radius, AR: $AR, N_eng: $N_eng, Vspec: $Vspec")
        push!(Vol_wing_list, 0.0)
        push!(Vol_nace_list, 0.0)

    end

    # error("message")  # NILS: stop here for debugging

end

ac_list_flat = vcat([ac_ref], ac_list)
colour_list = fill(colorant"#f66068", length(ac_list_flat))

p = PointLoadPlots.compare_stickfig(
    Vector{Any}(ac_list_flat); colour_list=colour_list, annotate_text = false, annotate_group = false, framestyle = :none,
)
# p = TASOPT.stickfig(ac_list_flat[4])
display(p)

# error("message")  # NILS: stop here for debugging

###############################################################################################

#                 end
#             end
#         end
#     end
# end
"""
error("message")  # NILS: stop here for debugging

### Save study results to .xlsx file for plotting in python

using XLSX
using DataFrames
using Tables

xlsx_filename = "volume_results_narrowbody_new.xlsx"

# Prepare vectors to collect flattened data
indep_var_flat = Vector{Any}()
CDS_flat = Vector{Any}()
WMTO_flat = Vector{Any}()
CDS_ref_flat = Vector{Any}()
WMTO_ref_flat = Vector{Any}()
Vol_wing_flat = Vector{Any}()

n1, n2, n3, n4, n5 = size(indep_var_grid)
for h in 1:n1, i in 1:n2, j in 1:n3, k in 1:n4, l in 1:n5
    indep = indep_var_grid[h, i, j, k, l]
    if isassigned(ac_grid, h, i, j, k, l)
        ac = ac_grid[h, i, j, k, l]
        if ac !== nothing && ac !== []
            push!(indep_var_flat, indep)
            push!(CDS_flat, ac.wing.layout.S * ac.para[iaCD, ipcruise1])
            push!(WMTO_flat, ac.parg[igWMTO])
            push!(CDS_ref_flat, ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1])
            push!(WMTO_ref_flat, ac_ref.parg[igWMTO])
            push!(Vol_wing_flat, 2 * (ac.wing.center.volume + ac.wing.inboard.volume + ac.wing.outboard.volume))
            continue
        end
    end
    # If not assigned or invalid, push NaN
    push!(indep_var_flat, indep)
    push!(CDS_flat, "NaN")
    push!(WMTO_flat, "NaN")
    push!(CDS_ref_flat, "NaN")
    push!(WMTO_ref_flat, "NaN")
    push!(Vol_wing_flat, "NaN")
end

# Split indep_var_flat into columns
fcs_loc_col = [x[1] for x in indep_var_flat]
radius_col = [x[2] for x in indep_var_flat]
AR_col = [x[3] for x in indep_var_flat]
N_eng_col = [x[4] for x in indep_var_flat]
Vspec_col = [x[5] for x in indep_var_flat]

df = DataFrame(
    fcs_loc = fcs_loc_col,
    radius = radius_col,
    AR = AR_col,
    N_eng = N_eng_col,
    Vspec = Vspec_col,
    CDS = CDS_flat,
    CDS_ref = CDS_ref_flat,
    WMTO = WMTO_flat,
    WMTO_ref = WMTO_ref_flat,
    Vol_wing = Vol_wing_flat,
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

# ac_ref = ac_grid[1, 1, 1, length(wing_frac_range)]
ac_1 = ac_grid[length(radius_range), 1, 1, 1]
ac_2 = ac_grid[1, 2, 1, 1]
ac_3 = ac_grid[1, 1, length(Neng_range), 1]
ac_4 = ac_grid[3, 3, 1, length(Vspec_range)]

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
# savefig(p, "stickfig_case_4.svg")
"""
