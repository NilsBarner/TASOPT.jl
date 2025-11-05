export update_fuse!, update_fuse_for_pax!
"""
update_fuse!(ac, imission)

Function to update the fuselage layout when there is a change in fuselage fuel-tank length.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac::aircraft`: object with aircraft parameters
    - `imission::Int64`: mission index

    **Outputs:**
    No direct outputs; parameters in `ac` are modified.
"""
function update_fuse!(ac, imission::Int64 = 1)
    #Unpack storage objects
    parg, parm, para, pare, options, fuse, fuse_tank, wing, htail, vtail, engine = unpack_ac(ac, imission) 
    parm = view(parm, :, imission)
    para = view(para, :, :, imission)

    nftanks = fuse_tank.tank_count #number of tanks
    lftank = parg[iglftank] # Get length of tank from previous iteration
    lftoffset = 2.0*ft_to_m #1 ft buffer for front and back of tanks

    #Useful relative distances to conserve
    lcyl = fuse.layout.l_cabin_cylinder
    dxeng2wbox = parg[igdxeng2wbox]
    dxapu2end = fuse.layout.x_end - fuse.APU.x
    dxshell2conend =fuse.layout.x_cone_end - fuse.layout.x_pressure_shell_aft
    dxshell2apu = fuse.APU.x - fuse.layout.x_pressure_shell_aft
    dxhbox2conend = fuse.layout.x_cone_end - htail.layout.box_x
    dxvbox2conend = fuse.layout.x_cone_end - vtail.layout.box_x

    if parg[igxftankaft] == 0.0 #if there is not a rear tank
        dxcyl2shellaft = fuse.layout.x_pressure_shell_aft - fuse.layout.x_end_cylinder
    else #if there is a rear tank
        dxcyl2shellaft = 0.0 #no need for offset between shell2 and blend2 since rear space cannot be used
    end

    #Update positions and fuselage length
    fuse.layout.x_end_cylinder = fuse.layout.x_start_cylinder + nftanks * (lftank + lftoffset) + lcyl

    fuse.layout.x_pressure_shell_aft = fuse.layout.x_end_cylinder + dxcyl2shellaft

    fuse.layout.x_cone_end = fuse.layout.x_pressure_shell_aft + dxshell2conend
    fuse.APU.r = [fuse.layout.x_pressure_shell_aft + dxshell2apu, 0.0, 0.0]
    fuse.layout.x_end = fuse.APU.x + dxapu2end
    fuse.HPE_sys.r = [fuse.layout.x_cone_end * 0.52484, 0.0,0.0]#TODO: address this
    
    htail.layout.box_x = fuse.layout.x_cone_end - dxhbox2conend
    vtail.layout.box_x = fuse.layout.x_cone_end - dxvbox2conend
    
    parg[igxeng    ] =  wing.layout.box_x - dxeng2wbox


    #Update fuselage aerodynamic parameters
    fuselage_drag!(fuse, parm, para, ipcruise1) #Recalculate fuselage bl properties

    #Update fuselage BL properties
    # Kinetic energy area at T.E.
    KAfTE = para[iaKAfTE, ipcruise1]
    # Surface dissapation area 
    DAfsurf = para[iaDAfsurf, ipcruise1]
    # Wake dissapation area
    DAfwake = para[iaDAfwake, ipcruise1]
    # Momentum area at ∞
    PAfinf = para[iaPAfinf, ipcruise1]

    # Assume K.E., Disspation and momentum areas are const. for all mission points:
    para[iaKAfTE, :] .= KAfTE
    para[iaDAfsurf, :] .= DAfsurf
    para[iaDAfwake, :] .= DAfwake
    para[iaPAfinf, :] .= PAfinf
end

"""
    update_fuse_for_pax!(ac)

Function to update the fuselage layout when the cabin length is not known a priori, for example if the radius is changed. 
It sizes the cabin for the design number of passengers.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `ac`::aircraft: aircraft object

    **Outputs:**
    Parameters in `parg` are modified. It also outputs:
    - `seats_per_row::Float64`: number of seats per row in main cabin (lower deck if double decker)
"""
function update_fuse_for_pax!(
    ac; run_from_within::Bool = false, lcyl::Float64 = 0.0, A_cargo::Float64 = 0.0,  # NILS: note the semicolon separating positional from keyword arguments
)  # last three arguments added by NILS
    parg, options, fuse, fuse_tank, wing, htail, vtail, _ = unpack_ac_components(ac)

    seat_pitch = fuse.cabin.seat_pitch
    seat_width = fuse.cabin.seat_width 
    aisle_halfwidth = fuse.cabin.aisle_halfwidth
    h_seat = fuse.cabin.seat_height
    d_floor = fuse.cabin.floor_distance
    front_seat_offset = fuse.cabin.front_seat_offset

    Rfuse = fuse.layout.radius
    dRfuse = fuse.layout.bubble_lower_downward_shift
    wfb = fuse.layout.bubble_center_y_offset
    nfweb = fuse.layout.n_webs

    # Calculate length occupied by FCS
    # Vfcs = parg[igVfcsfus]
    fcs_loc = parg[igfcsfusloc]
    Afcs = 0.0
    Acargo = 0.0
    if fcs_loc == 1.0  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl"
        lfcs = parg[iglfcsfus]  # NILS
        Vfcs = lfcs * pi * Rfuse^2  # NILS
    elseif fcs_loc == 0.0  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl
        lfcs = 0.0
        θ = -parg[igthetafloor]  # NILS (NOTE THE "-" SIGN)
    end

    #Find cabin length by placing seats
    if fuse.n_decks == 2 #if the aircraft is a double decker
        xopt, seats_per_row = optimize_double_decker_cabin(fuse) #Optimize the floor layout and passenger distributions

        lcyl, _ = find_double_decker_cabin_length(xopt, fuse) #Total length is maximum of the two

        paxsize = xopt[1]
        #Store angles
        fuse.cabin.floor_angle_main = xopt[2]
        fuse.cabin.floor_angle_top = find_floor_angles(true, Rfuse, dRfuse, θ1 = xopt[2], h_seat = h_seat, d_floor = d_floor)[2]
    else
        # NILS: find change in floor angle to accommodate underfloor FCS
        if fcs_loc == 1.0  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl
            θ = find_floor_angles(false, Rfuse, dRfuse, h_seat = h_seat) #Find the floor angle
            pim2θ = pi - 2 * (-θ)  # NILS (NOTE the "-" sign - look inside theta_from_area())
            A_cargo = Rfuse^2 / 2 * (pim2θ - sin(pim2θ))  # NILS
        elseif fcs_loc == 0.0  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl
            # NILS: calculate available cross-sectional area for FCS and cargo given θ
            ULD = fuse.cabin.unit_load_device
            ULDdims = UnitLoadDeviceDimensions[ULD]
            minwidth = ULDdims[2] #Base width
            c_ULD = minwidth
            h_ULD = ULDdims[1]
            dh_ULD = Rfuse - sqrt(Rfuse^2 - c_ULD^2 / 4)
            θfcs = -asin((Rfuse - dh_ULD - h_ULD) / Rfuse)
            Afcs = Rfuse^2 / 2 * ((pi + 2 * θ) - sin(pi + 2 * θ)) - (Rfuse^2 / 2 * ((pi + 2 * θfcs) - sin(pi + 2 * θfcs)))
            Acargo = Rfuse^2 / 2 * ((pi + 2 * θfcs) - sin(pi + 2 * θfcs))
            # println("c_ULD, h_ULD, dh_ULD, θfcs, Afcs, Acargo = $c_ULD, $h_ULD, $dh_ULD, $θfcs, $Afcs, $Acargo")
        end
        paxsize = fuse.cabin.exit_limit #maximum number of passengers
        w = find_cabin_width(Rfuse, wfb, nfweb, θ, h_seat) #Cabin width
        lcyl, _, seats_per_row = place_cabin_seats(paxsize, w, seat_pitch = seat_pitch, seat_width = seat_width, 
            aisle_halfwidth = aisle_halfwidth, front_seat_offset = front_seat_offset) #Cabin length
        # NILS: calculate volume from cross-sectional areas and cylinder length
        if fcs_loc == 0.0
            Vfcs = Afcs * lcyl
            Vcargo = Acargo * lcyl
        end
    end

    # Log available FCS volume
    parg[igVfcsfus] = Vfcs  # NILS

    #Useful relative distances to conserve
    dxeng2wbox = parg[igdxeng2wbox] #Distance from engine to wingbox
    dxcyl2shellaft = fuse.layout.x_pressure_shell_aft - fuse.layout.x_end_cylinder #Distance from blend2 to shell2
    dxapu2end = fuse.layout.x_end - fuse.APU.x #Distance from APU to end
    dxshell2conend = fuse.layout.x_cone_end - fuse.layout.x_pressure_shell_aft #Distance from shell2 to conend
    dxshell2apu = fuse.APU.x - fuse.layout.x_pressure_shell_aft #Distance from shell2 to APU
    dxhbox2conend = fuse.layout.x_cone_end - htail.layout.box_x #Distance from conend to xhbox
    dxvbox2conend = fuse.layout.x_cone_end - vtail.layout.box_x #Distance from conend to xvbox
    #Fraction of cabin length at which wing is located
    wbox_cabin_frac =  (wing.layout.box_x- fuse.layout.x_start_cylinder )/(fuse.layout.x_end_cylinder - fuse.layout.x_start_cylinder) 

    #######################################################################
    #Update positions and fuselage length
    # lcyl = lcyl + fuse.cabin.rear_seat_offset #Make cabin longer to leave room in the back
    lcyl = lcyl + fuse.cabin.rear_seat_offset + lfcs #Make cabin longer to leave room in the back  # lfcs added by NILS
    #######################################################################

    #Update positions and fuselage length
    #######################################################################
    # fuse.layout.x_end_cylinder = fuse.layout.x_start_cylinder + lcyl + lfcs  # lfcs added by NILS
    fuse.layout.x_end_cylinder = fuse.layout.x_start_cylinder + lcyl
    if fcs_loc == 1.0
        fuse.layout.x_centroid_fcs = fuse.layout.x_end_cylinder - lfcs / 2  # NILS: added to log longitudinal location of FCS centroid
    elseif fcs_loc == 0.0
        fuse.layout.x_centroid_fcs = fuse.layout.x_start_cylinder + lcyl / 2  # NILS: added to log longitudinal location of FCS centroid
    end
    #######################################################################

    #Update wingbox position
    wing.layout.box_x = fuse.layout.x_start_cylinder + wbox_cabin_frac * lcyl
       
    #Update other lengths
    fuse.layout.x_pressure_shell_aft = fuse.layout.x_end_cylinder + dxcyl2shellaft

    fuse.layout.x_cone_end = fuse.layout.x_pressure_shell_aft + dxshell2conend
    fuse.APU.r = [fuse.layout.x_pressure_shell_aft + dxshell2apu, 0.0,0.0]
    #######################################################################
    fuse.layout.x_end = fuse.APU.x + dxapu2end# + lfcs  # lfcs added by NILS
    #######################################################################
    fuse.HPE_sys.r = [fuse.layout.x_cone_end * 0.52484, 0.0, 0.0] #TODO: address this
    
    htail.layout.box_x = fuse.layout.x_cone_end - dxhbox2conend
    vtail.layout.box_x = fuse.layout.x_cone_end - dxvbox2conend
    
    parg[igxeng    ] =  wing.layout.box_x - dxeng2wbox #Move engine

    fuse.layout.l_cabin_cylinder = lcyl #Store new cabin length

    EvaluateCabinProps!(fuse) #Update cabin parameters

    # println("wcabin, θ = $wcabin, $θ")
    # println("lfcs = $lfcs")
    # println("lcyl = $lcyl")
    # println("seats_per_row = $seats_per_row")
    # println("fuse.layout.radius = $(fuse.layout.radius)")

    # # NILS: re-run function
    # if (Vfcs > 0.0 && fcs_loc == 0.0 && run_from_within == true)  # NILS: fcs_loc == 0.0 stands for "underfloor", whereas 1.0 stands for "rear" to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl
    #     println("Hello world")
    #     seats_per_row = update_fuse_for_pax!(
    #         ac; run_from_within=run_from_within, lcyl=lcyl, A_cargo=A_cargo,  # NILS: note the semicolon separating positional from keyword arguments
    #     )
    # end

    return seats_per_row
end

"""
    find_minimum_radius_for_seats_per_row(seats_per_row, ac_base)

This function calculates the minimum radius required to have a desired number of seats per row in the main cabin.

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `seats_per_row::Int64`: number of seats per row in main cabin (lower deck if double decker)
    - `ac_base::aircraft`: aircraft object

    **Outputs:**
    - `R::Float64`: minimum radius for desired number of seats per row (m)
"""
function find_minimum_radius_for_seats_per_row(seats_per_row::Int64, ac_base)
    ac = deepcopy(ac_base) #Copy input ac to avoid modifying it
    obj(x, grad) = x[1] + 1e3 * abs(check_seats_per_row_diff(seats_per_row, x, ac))  #Objective function is the radius plus a big penalty if constraint is not met

    #First, use global optimizer to find region of global optimum
    initial_x = [4.0]
    opt = Opt(:GN_DIRECT, length(initial_x)) #Use a global optimizer
    opt.lower_bounds = [0.0]
    opt.upper_bounds = [5.0]

    # opt_local = Opt(:GN_DIRECT, length(initial_x))
    opt.maxeval = 500  # Set the max number of evaluations
    # opt.local_optimizer = opt_local

    opt.min_objective = obj

    #Apply the equality constraint that ensures that the number of seats per row is the desired one
    #equality_constraint!(opt, (x, grad) -> check_seats_per_row_diff(seats_per_row, x, ac), 1e-5) 
    
    (minf,xopt,ret) = NLopt.optimize(opt, initial_x) #Solve optimization problem

    #Next, use local optimizer to polish off optimum
    opt = Opt(:LN_NELDERMEAD, length(initial_x)) #Use a local optimizer
    opt.lower_bounds = [0.0]
    opt.upper_bounds = [5.0]
    opt.min_objective = obj
    opt.ftol_rel = 1e-4

    (minf,xopt,ret) = NLopt.optimize(opt, xopt) #Solve optimization problem starting from global solution

    R = xopt[1]
    #Check if constraint is met
    diff = check_seats_per_row_diff(seats_per_row, xopt, ac)
    if diff ≈ 0.0
        return R
    else
        error("Optimizer failed to find a fuselage radius for the desired pax per row")
    end
end

"""
    check_seats_per_row_diff(seats_per_row, x, ac)

This function returns the difference between the desired number of seats per row and the one corresponding to a 
given radius

!!! details "🔃 Inputs and Outputs"
    **Inputs:**
    - `seats_per_row::Float64`: desired number of seats per row in main cabin (lower deck if double decker)
    - `x::Vector{Float64}`: vector with one entry containing the fuselage radius (m)
    - `ac::aircraft`: aircraft object

    **Outputs:**
    - `diff::Float64`: difference between desired number of seats per row and that for the input radius
"""
function check_seats_per_row_diff(seats_per_row, x, ac)
    Rfuse = x[1]
    ac.fuselage.layout.cross_section.radius = Rfuse
    try #Sometimes update_fuse_for_pax may fail
        seats_per_row_rad = update_fuse_for_pax!(ac)
        diff = seats_per_row_rad - seats_per_row
        #println("R = $Rfuse, s = $seats_per_row_rad")
        return diff
    catch
        return 1.0
    end
end