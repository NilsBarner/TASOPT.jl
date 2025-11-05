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
N_points = 20
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

# fcs_locs = ["nacelle", "wing", "fuselage"]
fcs_locs = ["wing"]
radius_range = range(80, 130, length=N_points) * 0.0254  # [m]
AR_range = range(6, 18, length=N_points)
tdivc_scale_range = range(0.75, 1.25, length=N_points)
Neng_range = [2, 4, 6, 8]
HTRf_range = range(0.2, 0.9, length=N_points)
# Vspec_range = range(3, 5.5, length=N_points) * 1e3 * 1e3  # [m^3]  Table 2 in FZO-PPN-REP-0031 - Fuel Cells Technical Report (STACK)
# Vspec_range = range(0.2, 0.4, length=N_points) * 1e3 * 1e3  # [m^3]  Table 2 in FZO-PPN-REP-0031 - Fuel Cells Technical Report (SYSTEM)
Vspec_range = range(0.2, 1.0, length=N_points) * 1e3 * 1e3  # NILS' estimate
fcs_fuselage_locations = [0.0, 1.0]

# # Investigate fuselage only
# fcs_locs = ["fuselage"]
# radius_range = range(80, 130, length=N_points) * 0.0254  # [m]
# AR_range = [10.1]
# tdivc_scale_range = [1.0]
# Neng_range = [2]
# HTRf_range = [0.3]
# Vspec_range = range(0.2, 1.0, length=N_points) * 1e3 * 1e3  # NILS' estimate
# fcs_fuselage_locations = [0.0, 1.0]

# # Investigate wing only
# fcs_locs = ["wing"]
# radius_range = 100 * 0.0254  # [m]
# AR_range = range(6, 18, length=N_points)
# tdivc_scale_range = range(0.75, 1.5, length=N_points)
# Neng_range = [2]
# HTRf_range = [0.3]
# Vspec_range = [0.5e6]
# fcs_fuselage_locations = [0.0]

# # Investigate nacelles only
# fcs_locs = ["nacelle"]
# radius_range = 100 * 0.0254  # [m]
# AR_range = [10.1]
# tdivc_scale_range = [1.0]
# Neng_range = [2, 4, 6, 8]
# HTRf_range = range(0.2, 0.9, length=N_points)
# Vspec_range = [0.5e6]
# fcs_fuselage_locations = [0.0]

indep_var_grid = Array{Any}(
    undef, length(fcs_locs), length(radius_range), length(AR_range), length(tdivc_scale_range), length(Neng_range), length(HTRf_range), length(Vspec_range), length(fcs_fuselage_locations)
)
ac_grid = Array{Any}(
    undef, length(fcs_locs), length(radius_range), length(AR_range), length(tdivc_scale_range), length(Neng_range), length(HTRf_range), length(Vspec_range), length(fcs_fuselage_locations)
)

for (h, fcs_loc) in enumerate(fcs_locs)
    for (i, radius) in enumerate(radius_range)
        for (j, AR) in enumerate(AR_range)
            for (k, tdivc_scale) in enumerate(tdivc_scale_range)
                for (l, N_eng) in enumerate(Neng_range)
                    for (m, HTR_f) in enumerate(HTRf_range)
                        for (n, Vspec) in enumerate(Vspec_range)
                            for (o, fcs_fuselage_location) in enumerate(fcs_fuselage_locations)

                                if (h == 1 && i == 1 && j == 1 && k == 1 && n == 1 && o == 1)
                                    radius = 100 * 0.0254  # [m]
                                    AR = 10.1
                                    tdivc_scale = 1.0
                                    Vspec = 0.5e6
                                    fcs_fuselage_location = 0.0
                                elseif (h == 2 && i == 1 && l == 1 && m == 1 && n == 1 && o == 1)
                                    radius = 100 * 0.0254  # [m]
                                    N_eng = 2
                                    HTR_f = 0.3
                                    Vspec = 0.5e6
                                    fcs_fuselage_location = 0.0
                                elseif (h == 3 && j == 1 && k == 1 && l == 1 && m == 1)
                                    AR = 10.1
                                    tdivc_scale = 1.0
                                    N_eng = 2
                                    HTR_f = 0.3
                                end

                                indep_var_grid[h, i, j, k, l, m, n, o] = (fcs_loc, radius, AR, tdivc_scale, N_eng, HTR_f, Vspec, fcs_fuselage_location)

                                # if fcs_loc == "fuselage" || fcs_fuselage_location != 1.0

                                if (h == 1 && i == 1 && j == 1 && k == 1 && n == 1 && o == 1) ||
                                   (h == 2 && i == 1 && l == 1 && m == 1 && n == 1 && o == 1) ||
                                   (h == 3 && j == 1 && k == 1 && l == 1 && m == 1)

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
                                        size_aircraft!(ac, iter=50, printiter=false)  # disabled default iteration print-outs to console
                                        # Store aircraft object in output grids
                                        ac_grid[h, i, j, k, l, m, n, o] = ac

                                    catch e
                                        println("fcs_loc: $fcs_loc, radius: $radius, AR: $AR, tdivc_scale: $tdivc_scale, N_eng: $N_eng, HTR_f: $HTR_f, Vspec: $Vspec, fcs_fuselage_location: $fcs_fuselage_location")
                                    end
                                    
                                end

                                # error("message")  # NILS: stop here for debugging
                            end
                        end
                    end
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

xlsx_filename = "volume_results_narrowbody_wing_new.xlsx"

# Prepare vectors to collect flattened data
indices_flat = Vector{Any}()
indep_var_flat = Vector{Any}()
CDS_flat = Vector{Any}()
WMTO_flat = Vector{Any}()
CDS_ref_flat = Vector{Any}()
WMTO_ref_flat = Vector{Any}()
Vol_wing_flat = Vector{Any}()
PFEI_flat = Vector{Any}()
PFEI_ref_flat = Vector{Any}()
seats_abreast_flat = Vector{Any}()
seats_abreast_ref_flat = Vector{Any}()
Vol_nacelle_flat = Vector{Any}()

n1, n2, n3, n4, n5, n6, n7, n8 = size(indep_var_grid)
for h in 1:n1, i in 1:n2, j in 1:n3, k in 1:n4, l in 1:n5, m in 1:n6, n in 1:n7, o in 1:n8
    indep = indep_var_grid[h, i, j, k, l, m, n, o]
    if isassigned(ac_grid, h, i, j, k, l, m, n, o)
        ac = ac_grid[h, i, j, k, l, m, n, o]
        if ac !== nothing && ac !== []
            push!(indices_flat, (h, i, j, k, l, m, n, o))
            push!(indep_var_flat, indep)
            push!(CDS_flat, ac.wing.layout.S * ac.para[iaCD, ipcruise1])
            push!(WMTO_flat, ac.parg[igWMTO])
            push!(CDS_ref_flat, ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1])
            push!(WMTO_ref_flat, ac_ref.parg[igWMTO])
            push!(Vol_wing_flat, 2 * (ac.wing.center.volume + ac.wing.inboard.volume + ac.wing.outboard.volume))
            push!(PFEI_flat, ac.parm[imPFEI])
            push!(PFEI_ref_flat, ac_ref.parm[imPFEI])
            push!(seats_abreast_flat, ac.fuselage.cabin.seats_abreast_main)
            push!(seats_abreast_ref_flat, ac_ref.fuselage.cabin.seats_abreast_main)
            push!(Vol_nacelle_flat, ac.parg[igVfcsavnacetot])
            continue
        end
    end
    # # If not assigned or invalid, push NaN
    # push!(indices_flat, (h, i, j, k, l, m, n, o))
    # push!(indep_var_flat, indep)
    # push!(CDS_flat, "NaN")
    # push!(WMTO_flat, "NaN")
    # push!(CDS_ref_flat, "NaN")
    # push!(WMTO_ref_flat, "NaN")
    # push!(Vol_wing_flat, "NaN")
    # push!(PFEI_flat, "NaN")
    # push!(PFEI_ref_flat, "NaN")
    # push!(seats_abreast_flat, "NaN")
    # push!(seats_abreast_ref_flat, "NaN")
    # push!(Vol_nacelle_flat, "NaN")
end

# Split indep_var_flat into columns
index_col = [join(string.(inds), "_") for inds in indices_flat]
fcs_loc_col = [x[1] for x in indep_var_flat]
radius_col = [x[2] for x in indep_var_flat]
AR_col = [x[3] for x in indep_var_flat]
tdivc_scale_col = [x[4] for x in indep_var_flat]
N_eng_col = [x[5] for x in indep_var_flat]
HTR_f_col = [x[6] for x in indep_var_flat]
Vspec_col = [x[7] for x in indep_var_flat]
fcs_fuselage_location_col = [x[8] for x in indep_var_flat]

df = DataFrame(
    index = index_col,
    fcs_loc = fcs_loc_col,
    radius = radius_col,
    AR = AR_col,
    tdivc_scale = tdivc_scale_col,
    N_eng = N_eng_col,
    HTR_f = HTR_f_col,
    Vspec = Vspec_col,
    fcs_fuselage_location = fcs_fuselage_location_col,
    CDS = CDS_flat,
    CDS_ref = CDS_ref_flat,
    WMTO = WMTO_flat,
    WMTO_ref = WMTO_ref_flat,
    Vol_wing = Vol_wing_flat,
    PFEI = PFEI_flat,
    PFEI_ref = PFEI_ref_flat,
    seats_abreast = seats_abreast_flat,
    seats_abreast_ref = seats_abreast_ref_flat,
    Vol_nacelle = Vol_nacelle_flat,
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

