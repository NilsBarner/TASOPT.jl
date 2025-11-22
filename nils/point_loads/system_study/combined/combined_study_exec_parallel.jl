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
include(joinpath(pwd(), "nils/point_loads/system_study/volume/bin_packing/coordinate_compilation/wingbox_coordinate_compiler.jl"))

# ---------------- Parallel-run support ----------------
# Accept an optional sigma value (single) and optional suffix on the command line.
# Usage:
#   julia --project=. script.jl 1500 mysuffix
# or to run the whole sweep (old behaviour):
#   julia --project=. script.jl

sigma_arg = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : nothing
out_suffix = length(ARGS) >= 2 ? ARGS[2] : (sigma_arg !== nothing ? string(Int(round(sigma_arg))) : "allsigs")

if sigma_arg !== nothing
    sigma_fcs_range = [sigma_arg]
else
    sigma_fcs_range = collect(range(1.5, 4, length=10)) * 1e3
end

# canonical toml path for this job (unique per suffix)
function job_toml_path(ac_segment::String, suffix::AbstractString)
    if ac_segment == "narrowbody"
        return joinpath(pwd(), "example", "a220100_lh2_input_$(suffix).toml")
    elseif ac_segment == "regional"
        return joinpath(pwd(), "example", "atr72600_lh2_input_$(suffix).toml")
    end
end
# ---------------- end parallel-run support ----------------

###

# Top-level settings
ac_segment = "narrowbody"  # "narrowbody" or "regional"
##### NILS: RESET "REGIONAL" - CURRENTLY DOES NOT RUN, PRESUMABLY BECAUSE OF INFEASIBLE INPUT FROM MASS-EXPERIMENTS #####
N_points = 19

# ##############################

# datafile = joinpath(pwd(), "example", "cryo_input.toml")
# templatefile = joinpath(TASOPT.__TASOPTroot__, "../example/defaults/default_input.toml")
# data = TOML.parsefile(datafile)
# default = TOML.parsefile(templatefile)
# fuse = read_input("Fuselage", data, default)
# dfuse = default["Fuselage"]
# geom = read_input("Geometry", fuse, dfuse)
# dgeom = dfuse["Geometry"]
# readgeom(x) = read_input(x, geom, dgeom)
# x_start_cylinder = Distance(readgeom("x_start_cylinder"))
# println("x_start_cylinder")
# error("stop")
# ##############################

###

# parse index column if it was read from CSV as a string like "(1, 9, 8, 10, 1, 1, 1, 7, 2, 1)"
function parse_ntuple10(s::AbstractString)
    s2 = strip(s)
    s2 = s2[2:end-1]
    parts = split(s2, ",")
    vals = parse.(Int64, strip.(parts))
    return convert(NTuple{10,Int64}, tuple(vals...))
end

###

function analyse_reference_aircraft(ac_segment::String, toml_path::AbstractString)

    # Baseline LH2 aircraft

    # Set all point loads applied via .toml file to 0
    toml_data = TOML.parsefile(toml_path)   # parse the job-specific toml (if missing, TOML.parsefile will fail; next lines overwrite anyway)
    # toml_data["Propulsion"]["number_of_engines"] = 2
    # toml_data["Propulsion"]["Turbomachinery"]["HTR_fan"] = 0.3
    # toml_data["Fuselage"]["l_fcs_fuselage"] = 0.0
    # toml_data["Fuselage"]["Geometry"]["unit_load_device"] = "LD3-45"
    # toml_data["Fuselage"]["Geometry"]["theta_floor"] = 0.0
    # toml_data["Fuselage"]["fcs_fuselage_location"] = 1.0
    # toml_data["Wing"]["AR"] = 10.1
    # toml_data["Wing"]["tdivc_scale"] = 1.0
    # toml_data["Wing"]["TR_scale"] = 1.0
    # toml_data["Wing"]["has_strut"] = false
    # toml_data["Fuselage"]["Geometry"]["radius"] = 100 * 0.0254
    # toml_data["Propulsion"]["Weight"]["custom_weight_delta"] = 0.0
    # # write job-specific toml
    # open(toml_path, "w") do file
    #     TOML.print(file, toml_data)
    # end

    # Read model inputs from that same job-specific toml
    global ac_ref = read_aircraft_model(toml_path)

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
    x_centroid_fuse = ac_ref.fuselage.layout.x_centroid_fcs

    # Calculate maximum propulsive power throughout mission
    T_net_max_ref = ac_ref.pare[ieFe, :] * ac_ref.parg[igneng]  # NILS: see calculate_thrust_from_ROC!() for proof that pare[ieFe] stands for per-engine thrust
    V_0_max_ref = ac_ref.pare[ieu0, :]
    P_prop_max_ref = maximum(T_net_max_ref .* V_0_max_ref)

    df = DataFrame(
        index = 0,
        sigma_fcs = 0.0,
        fcs_loc = "NA",
        radius = toml_data["Fuselage"]["Geometry"]["radius"],
        AR = toml_data["Wing"]["AR"],
        tdivc_scale = toml_data["Wing"]["tdivc_scale"],
        N_eng = toml_data["Propulsion"]["number_of_engines"],
        HTR_f = toml_data["Propulsion"]["Turbomachinery"]["HTR_fan"],
        l_fcs = toml_data["Fuselage"]["l_fcs_fuselage"],
        theta_floor = toml_data["Fuselage"]["Geometry"]["theta_floor"],
        fcs_fuselage_location = toml_data["Fuselage"]["fcs_fuselage_location"],
        has_strut = toml_data["Wing"]["has_strut"],

        CDS = ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1],
        WMTO = ac_ref.parg[igWMTO],
        Vol_wing = 2*(ac_ref.wing.center.volume + ac_ref.wing.inboard.volume + ac_ref.wing.outboard.volume),
        PFEI = ac_ref.parm[imPFEI],
        seats_abreast = ac_ref.fuselage.cabin.seats_abreast_main,
        Vol_nacelle = ac_ref.nils.V_fcs_nacelle,  # parg[igVfcsavnacetot],
        y_centroid_wing = y_centroid_wing,
        L_fuse = ac_ref.fuselage.layout.l_cabin_cylinder,
        y_centroid_nacelles = y_centroid_nacelles,
        span = b,
        x_centroid_fuse = x_centroid_fuse,
        Vol_fuse = ac_ref.nils.V_fcs_fuselage,  # ac_ref.parg[igVfcsfus],
        length = ac_ref.fuselage.layout.x_end,
        P_prop_max = P_prop_max_ref,

        # Wto = ac_ref.parg[igWMTO],
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
        # CDS = ac_ref.wing.layout.S * ac_ref.para[iaCD, ipcruise1],
        CDfuseS = ac_ref.wing.layout.S * ac_ref.para[iaCDfuse, ipcruise1],
        CDiS = ac_ref.wing.layout.S * ac_ref.para[iaCDi, ipcruise1],
        CDwingS = ac_ref.wing.layout.S * ac_ref.para[iaCDwing, ipcruise1],
        CDhtailS = ac_ref.wing.layout.S * ac_ref.para[iaCDhtail, ipcruise1],
        CDvtailS = ac_ref.wing.layout.S * ac_ref.para[iaCDvtail, ipcruise1],
    )

    # println("ac_ref.pare = ", ac_ref.pare[ieFe, ipclimbn])
    # println("ac_ref.pare = ", ac_ref.pare[ieFe, :])
    # println("ac_ref.pare = ", ac_ref.pare[ieu0, :])

    # at end: write the reference csv with suffix to avoid collisions
    CSV.write("combined_results_narrowbody_ref_101125_$(out_suffix).csv", df)

    return ac_ref, Weng_single_ref

end

###

# function run_study()

# Model reference aircraft
toml_path = job_toml_path(ac_segment, out_suffix)
ac_ref, Weng_single_ref = analyse_reference_aircraft(ac_segment, toml_path)
println("Weng_single_ref =", Weng_single_ref)
# error("message")

###

# result vectors
index_vec = Vector{NTuple{10,Int}}()
sigma_fcs_vec = Float64[]
fcs_loc_vec = String[]
radius_vec = Float64[]
AR_vec = Float64[]
tdivc_vec = Float64[]
Neng_vec = Int[]
HTRf_vec = Float64[]
lfcs_vec = Float64[]
thetafloor_vec = Float64[]
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
x_centroid_fuse_vec = Float64[]
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
# CDS_vec = Float64[]
CDfuseS_vec = Float64[]
CDiS_vec = Float64[]
CDwingS_vec = Float64[]
CDhtailS_vec = Float64[]
CDvtailS_vec = Float64[]

#####

# df = CSV.read(joinpath(pwd(), "nils", "point_loads", "system_study", "volume", "volume_results_narrowbody_101125.csv"), DataFrame)
df = CSV.read(joinpath(pwd(), "volume_results_narrowbody_101125.csv"), DataFrame)
# P_fcs = 20e6
# sigma_fcs_range = collect(range(1.5, 4, length=10)) * 1e3
# # sigma_fcs_range = [2e3]
# sigma_fcs_range = [4e3]
# sigma_fcs_range = [1.5e3]

length_sigma_fcs_range = length(sigma_fcs_range)
length_df = length(eachrow(df))

# Data

sigma_gt = 8e3  # W/kg based on CFM Leap in GT_baseline_slides.pptx
sigma_fcs_wing_min = 2e3  # W/kg (entire FCS in nacelle)
sigma_fcs_nacelle_max = 12e3  # W/kg (only motor + inverter in nacelle)
sigma_fcs_min = 2e3  # W/kg (minimum for entire FCS in nacelle)
sigma_fcs_max = 4e3  # W/kg (maximum for entire FCS in nacelle)

counter_out = 0
for (i, sigma_fcs) in enumerate(sigma_fcs_range)
    
    counter_in = 0
    for (j, row) in enumerate(eachrow(df))
        
        # Convert index tuple from string to NTuple{10,Int}
        raw_index = row.index
        index = parse_ntuple10(raw_index)
        fcs_loc = row.fcs_loc
        radius = row.radius
        AR = row.AR
        tdivc_scale = row.tdivc_scale
        N_eng = row.N_eng
        HTR_f = row.HTR_f
        l_fcs = row.l_fcs
        theta_floor = row.theta_floor
        fcs_fuselage_location = row.fcs_fuselage_location
        has_strut = row.has_strut

        CDS = row.CDS
        WMTO = row.WMTO
        Vol_wing = row.Vol_wing
        PFEI = row.PFEI
        seats_abreast = row.seats_abreast
        Vol_nacelle = row.Vol_nacelle
        y_centroid_wing = row.y_centroid_wing
        L_fuse = row.L_fuse
        y_centroid_nacelles = row.y_centroid_nacelles
        span = row.span
        x_centroid_fuse = row.x_centroid_fuse
        Vol_fuse = row.Vol_fuse
        _length = row.length
        ###
        # root_span = ac.fuselage.root_span
        # dfan = ac.fuselage.dfan
        ###

        # Define wing and fuselage point loads

        W_fcs = Weng_single_ref * sigma_gt / sigma_fcs

        if fcs_loc == "nacelle"
            nacelle_frac = 1.0
        elseif (fcs_loc == "fuselage" || fcs_loc == "wing")
            nacelle_frac = 0.0
        end

        # sigma_fcs_wing = sigma_fcs_wing_max - nacelle_frac * (sigma_fcs_wing_max - sigma_fcs)
        # W_fcs_wing = Weng_single_ref * sigma_gt / sigma_fcs_wing
        # W_fcs_fus = W_fcs - W_fcs_wing

        sigma_fcs_nacelle = sigma_fcs_nacelle_max - nacelle_frac * (sigma_fcs_nacelle_max - sigma_fcs)
        W_fcs_nacelle = Weng_single_ref * sigma_gt / sigma_fcs_nacelle
        custom_weight_delta = W_fcs_nacelle - Weng_single_ref
        Fz_point_nacelle = -(W_fcs_nacelle - Weng_single_ref)

        ###

        # Write engine point loads to .toml file

        # use the job-specific toml path:
        toml_data = TOML.parsefile(toml_path)
        toml_data["Propulsion"]["number_of_engines"] = N_eng
        toml_data["Propulsion"]["Turbomachinery"]["HTR_fan"] = HTR_f
        if fcs_loc == "nacelle"
            toml_data["Fuselage"]["l_fcs_fuselage"] = 0.0
            toml_data["Fuselage"]["Geometry"]["theta_floor"] = 0.0  # [deg]
            toml_data["Wing"]["has_strut"] = false
            # b = span
            # bo = root_span
            # ηs = toml_data["Wing"]["panel_break_location"]
            # dy = 2 * ac.parg[igdfan]  # space to leave near wing root and tip [m]
            # if N_eng == 2
            #     yi = ηs * b / 2
            # else
            #     yi = LinRange(bo / 2 + dy, b / 2 * 3 / 4, Int(N_eng / 2))
            # end
        elseif fcs_loc == "wing"                                                
            toml_data["Fuselage"]["l_fcs_fuselage"] = 0.0
            toml_data["Fuselage"]["Geometry"]["theta_floor"] = 0.0  # [deg]
            toml_data["Wing"]["has_strut"] = has_strut
        elseif fcs_loc == "fuselage"
            toml_data["Fuselage"]["l_fcs_fuselage"] = l_fcs
            toml_data["Fuselage"]["Geometry"]["theta_floor"] = theta_floor  # [deg]
            toml_data["Fuselage"]["fcs_fuselage_location"] = fcs_fuselage_location  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl"
            toml_data["Wing"]["has_strut"] = false
        end
        toml_data["Wing"]["AR"] = AR
        toml_data["Wing"]["tdivc_scale"] = tdivc_scale
        toml_data["Wing"]["TR_scale"] = 1.0
        toml_data["Fuselage"]["Geometry"]["radius"] = radius

        # Apply point load in centroid of volume

        if fcs_loc == "nacelle"
            sigma_fcs_nacelle = sigma_fcs_nacelle_max - nacelle_frac * (sigma_fcs_nacelle_max - sigma_fcs)
            W_fcs_nacelle = Weng_single_ref * sigma_gt / sigma_fcs_nacelle
            custom_weight_delta = W_fcs_nacelle - Weng_single_ref
            println("custom_weight_delta =", -Fz_point_nacelle)
            Fz_point_nacelle = -(W_fcs_nacelle - Weng_single_ref)
            toml_data["Propulsion"]["Weight"]["custom_weight_delta"] = -Fz_point_nacelle  # point load at engine location on wing
        elseif fcs_loc == "wing"
            toml_data["Propulsion"]["Weight"]["custom_weight_delta"] = 0.0
        elseif fcs_loc == "fuselage"
            toml_data["Propulsion"]["Weight"]["custom_weight_delta"] = 0.0
        end

        open(toml_path, "w") do file
            TOML.print(file, toml_data)
        end

        # read the model using the same toml_path
        global ac = read_aircraft_model(toml_path)

        if fcs_loc == "nacelle"
            nothing

        elseif fcs_loc == "wing"

            frac_span = y_centroid_wing / (span / 2)
            span_loc = Dict("frac_span" => frac_span)
            Fz_point_wing = -W_fcs
            # println("span_loc, Fz_point_wing =", span_loc, ", ", Fz_point_wing)

            # Add wing point loads to mutable structure
            wing_load = TASOPT.structures.PointLoad(
                force = SVector(0.0, 0.0, Fz_point_wing),
                r = SVector(0.0, span_loc, 0.0),
                frame = TASOPT.WORLD
            )
            TASOPT.structures.add_wing_point_load!(ac.wing, wing_load)

            # Add fuselage point loads to mutable structure
            fus_load = TASOPT.structures.PointLoad(
                force = SVector(0.0, 0.0, 0.0),  # NILS: previously, I mistakenly multiplied by number of engines here (20.10.2025)
                r = SVector(0.0, 0.0, 0.0),
                frame = TASOPT.WORLD
            )
            TASOPT.structures.add_fus_point_load!(ac.fuselage, fus_load)

        elseif fcs_loc == "fuselage"

            x_start_cylinder = 20 * 0.3048  # TO-DO @NILS: read from toml_data
            # x_start_cylinder = toml_data["Fuselage"]["Geometry"]["x_start_cylinder"]
            frac_fuse = (x_centroid_fuse - x_start_cylinder) / L_fuse
            _fcs_loc = Dict("frac_len" => frac_fuse)
            Fz_point_fus = -W_fcs
            
            # Add fuselage point loads to mutable structure
            fus_load = TASOPT.structures.PointLoad(
                force = SVector(0.0, 0.0, Fz_point_fus),  # NILS: previously, I mistakenly multiplied by number of engines here (20.10.2025)
                r = SVector(_fcs_loc, 0.0, 0.0),
                frame = TASOPT.WORLD
            )
            TASOPT.structures.add_fus_point_load!(ac.fuselage, fus_load)

            # Add wing point loads to mutable structure
            wing_load = TASOPT.structures.PointLoad(
                force = SVector(0.0, 0.0, 0.0),
                r = SVector(0.0, 0.0, 0.0),
                frame = TASOPT.WORLD
            )
            TASOPT.structures.add_wing_point_load!(ac.wing, wing_load)

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
            x_centroid_fuse = ac.fuselage.layout.x_centroid_fcs

            # Calculate maximum propulsive power throughout mission
            T_net_max = ac.pare[ieFe, :] * ac.parg[igneng]  # NILS: see calculate_thrust_from_ROC!() for proof that pare[ieFe] stands for per-engine thrust
            V_0_max = ac.pare[ieu0, :]
            P_prop_max = maximum(T_net_max .* V_0_max)
            
            # save results immediately
            push!(index_vec, index)
            push!(sigma_fcs_vec, sigma_fcs)
            push!(fcs_loc_vec, fcs_loc)
            println("sigma_fcs, fcs_loc =", sigma_fcs, ", ", fcs_loc)
            push!(radius_vec, radius)
            push!(AR_vec, AR)
            push!(tdivc_vec, tdivc_scale)
            push!(Neng_vec, N_eng)
            push!(HTRf_vec, HTR_f)
            push!(lfcs_vec, l_fcs)
            push!(thetafloor_vec, theta_floor)
            push!(fcs_fuselage_loc_vec, fcs_fuselage_location)
            push!(strut_bool_vec, has_strut)

            push!(CDS_vec, ac.wing.layout.S * ac.para[iaCD, ipcruise1])
            push!(WMTO_vec, ac.parg[igWMTO])
            push!(Vol_wing_vec, 2*(ac.wing.center.volume + ac.wing.inboard.volume + ac.wing.outboard.volume))
            push!(PFEI_vec, ac.parm[imPFEI])
            push!(seats_abreast_vec, ac.fuselage.cabin.seats_abreast_main)
            push!(Vol_nacelle_vec, ac.nils.V_fcs_nacelle)  # ac.parg[igVfcsavnacetot])
            push!(y_centroid_wing_vec, y_centroid_wing)
            push!(L_fuse_vec, ac.fuselage.layout.l_cabin_cylinder)
            push!(y_centroid_nacelles_vec, y_centroid_nacelles)
            push!(span_vec, b)
            push!(x_centroid_fuse_vec, x_centroid_fuse)
            push!(Vol_fuse_vec, ac.nils.V_fcs_fuselage)  # ac.parg[igVfcsfus])
            push!(length_vec, ac.fuselage.layout.x_end)
            push!(P_prop_max_vec, P_prop_max)

            # push!(Wto_vec, ac.parg[igWMTO])
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
            # push!(CDS_vec, ac.wing.layout.S * ac.para[iaCD, ipcruise1])
            push!(CDfuseS_vec, ac.wing.layout.S * ac.para[iaCDfuse, ipcruise1])
            push!(CDiS_vec, ac.wing.layout.S * ac.para[iaCDi, ipcruise1])
            push!(CDwingS_vec, ac.wing.layout.S * ac.para[iaCDwing, ipcruise1])
            push!(CDhtailS_vec, ac.wing.layout.S * ac.para[iaCDhtail, ipcruise1])
            push!(CDvtailS_vec, ac.wing.layout.S * ac.para[iaCDvtail, ipcruise1])

        catch e
            # @warn "failed case" fcs_loc=fcs_loc radius=radius AR=AR exception=e
            # println("fcs_loc: $fcs_loc, radius: $radius, AR: $AR, tdivc_scale: $tdivc_scale, N_eng: $N_eng, HTR_f: $HTR_f, l_fcs: $l_fcs, theta_floor: $theta_floor, fcs_fuselage_location: $fcs_fuselage_location", " has_strut: $has_strut")
            nothing
        end

        counter_in += 1
        println("Inner progress: ", counter_in, " / ", length_df, " | Outer progress: ", counter_out, " / ", length_sigma_fcs_range)
    end

    global counter_out += 1

end

df = DataFrame(
    index = index_vec,
    sigma_fcs = sigma_fcs_vec,
    fcs_loc = fcs_loc_vec,
    radius = radius_vec,
    AR = AR_vec,
    tdivc_scale = tdivc_vec,
    N_eng = Neng_vec,
    HTR_f = HTRf_vec,
    l_fcs = lfcs_vec,
    theta_floor = thetafloor_vec,
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
    x_centroid_fuse = x_centroid_fuse_vec,
    Vol_fuse = Vol_fuse_vec,
    length = length_vec,
    P_prop_max = P_prop_max_vec,

    # Wto = ac_ref.parg[igWMTO]
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
    # CDS = CDS_vec,
    CDfuseS = CDfuseS_vec,
    CDiS = CDiS_vec,
    CDwingS = CDwingS_vec,
    CDhtailS = CDhtailS_vec,
    CDvtailS = CDvtailS_vec,
)

# return df
# end

# df = run_study()

CSV.write("combined_results_narrowbody_101125_$(out_suffix).csv", df)
