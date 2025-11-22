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
    ac; Vfcs_is_input::Bool = false, run_from_within::Bool = false, lcyl::Float64 = 0.0, A_cargo::Float64 = 0.0,  # NILS: note the semicolon separating positional from keyword arguments
)  # last three arguments added by NILS
    """ NILS
    This function is called only once in read_input.jl prior to sizing the aircraft to
    establish the required cabin length `fuse.layout.l_cabin_cylinder` (`lcyl`)
    for a given payload-range requirement I prescribe, which is then kept constant,
    whereas update_fuse! is called multiple times during sizing to update the overall
    length of the cylindrical fuselage section
    `fuse.layout.x_end_cylinder = fuse.layout.x_start_cylinder + nftanks * (lftank + lftoffset) + lcyl`.

    NOTE 1: while it may seem stange that certain dimensional inputs like fuse.layout.x_end_cylinder or
    fuse.layout.x_end (see read_input.jl) are called in this function and never see the updated values during sizing,
    those are merely used to define relative distances to be conserved, like that between the APU
    and the fuselage end (`dxapu2end`). It is IMPORTANT though that the values of x_end_cylinder and
    x_end in any .toml files are sensible RELATIVE to each other, otherwise the aircraft geometry will be infeasible.

    NOTE 2: I only have to make updates to the fuselage geometry here, not in update_fuse!, because the FCS
    gets placed within the cabin section of the cylindrical fuselage (between fuse.layout.x_start_cylinder
    and the aft-fuselage tanks), which is carried into update_fuse! via `lcyl` (cabin length).
    """
    parg, options, fuse, fuse_tank, wing, htail, vtail, _ = unpack_ac_components(ac)

    seat_pitch = fuse.cabin.seat_pitch
    seat_width = fuse.cabin.seat_width 
    aisle_halfwidth = fuse.cabin.aisle_halfwidth
    h_seat = fuse.cabin.seat_height
    d_floor = fuse.cabin.floor_distance
    front_seat_offset = fuse.cabin.front_seat_offset

    Rfuse = fuse.layout.radius
    dRfuse = fuse.layout.bubble_lower_downward_shift  # NILS: zero in my case (see Figure 4 in tasopt.pdf)
    wfb = fuse.layout.bubble_center_y_offset  # NILS: zero in my case (see Figure 4 in tasopt.pdf)
    nfweb = fuse.layout.n_webs

    """ NILS
    Newton solver for pim2θ in pim2θ - sin pim2θ = K, K = 2S/r^2.
    In words, by how much must θ be increased to accommodate the
    FCS without shrinking the current cargo volume? For equal
    cargo and FCS compartment length, this translates into
    solving for the circular segment area Afcs + Acargo.
    # https://en.wikipedia.org/wiki/Circular_segment
    """
    function pim2θ_from_area(r, S; tol=1e-12, maxiter=50)
        S <= 0 && return 0.0  # trivial
        K = 2S/(r*r)
        pim2θ = (6K)^(1/3)  # very good initial guess for small K
        for i in 1:maxiter
            f = pim2θ - sin(pim2θ) - K
            abs(f) < tol && return pim2θ
            df = 1 - cos(pim2θ)
            if abs(df) < 1e-16
                pim2θ += 1e-6  # dodge near-zero derivative (rare)
            else
                pim2θ -= f/df
            end
        end
        return pim2θ
    end

    """ NILS
    Imports and initialisations. Vfcs_is_input == true
    is not being used currently (aim is to calculate
    required cabin floor angle and thus -length to
    maintain a constant cargo volume for a given FCS volume),
    but I am leaving this function argument in place for future use.

    fcs_loc == 0.0 stands for "underfloor", whereas
    fcs_loc == 1.0 stands for "rear"
    to comply with numeric slot `parg = zeros(Float64, igtotal)` in read_input.jl
    """
    fcs_loc = ac.nils.fcs_fuselage_location  # parg[igfcsfusloc]
    Afcs = 0.0
    Acargo = 0.0

    if Vfcs_is_input == true
        Vfcs = ac.nils.V_fcs_fuselage  # parg[igVfcsfus]
        if fcs_loc == 1.0
            lfcs = Vfcs / (pi * Rfuse^2)  # length of disk with radius Rfuse and volume Vfcs
        elseif fcs_loc == 0.0
            lfcs = 0.0  # included in cabin length
        end

    elseif Vfcs_is_input == false
        if fcs_loc == 1.0  
            lfcs = ac.nils.l_fcs_fuselage  # parg[iglfcsfus]
        elseif fcs_loc == 0.0
            lfcs = 0.0  # included in cabin length
            θ = -ac.nils.theta_floor  # parg[igthetafloor]
            """ NILS
            (NOTE the "-" to comply with the sign convention in
            https://mit-lae.github.io/TASOPT.jl/dev/structures/cabin_sizing/#Theory-Double-decker-aircraft-7a29cdcf253b565d ,
            where angles are measured positive upwards from the fuselage centreline)
            """
        end
    end

    #Find cabin length by placing seats
    if fuse.n_decks == 2 #if the aircraft is a double decker
        xopt, seats_per_row = optimize_double_decker_cabin(fuse) #Optimize the floor layout and passenger distributions

        lcyl, _ = find_double_decker_cabin_length(xopt, fuse) #Total length is maximum of the two

        paxsize = xopt[1]
        #Store angles
        fuse.cabin.floor_angle_main = xopt[2]
        fuse.cabin.floor_angle_top = find_floor_angles(true, Rfuse, dRfuse, θ1 = xopt[2], h_seat = h_seat, d_floor = d_floor)[2]

        """ NILS
        I moved the below code block from after this if-else statement to here,
        because for the single-decker case with no FCS it turned out to yield the
        same result as the first `find_floor_angles` in src/structures/update_fuse.jl
        in the original TASOPT.jl repo. Only if fuse.n_decks == 2
        then paxsize, which is an output of the optimisation, does not get used
        anywhere else but in the below code block, whereas elseif fuse.n_decks != 2
        paxsize already got computed, and d_floor, h_seat, bubble_lower_downward_shift,
        bubble_center_y_offset, and Rfuse are constants.
        """
        #Find new cabin length
        d_floor = fuse.cabin.floor_distance
        h_seat = fuse.cabin.seat_height 
        θ = find_floor_angles(false, fuse.layout.radius, fuse.layout.cross_section.bubble_lower_downward_shift, h_seat = h_seat, d_floor=d_floor) #Find the floor angle
        wcabin = find_cabin_width(fuse.layout.radius, fuse.layout.bubble_center_y_offset, fuse.layout.n_webs, θ, h_seat) #Find cabin width
        lcyl, _, _ = place_cabin_seats(paxsize, wcabin, seat_pitch = seat_pitch, seat_width = seat_width, 
            aisle_halfwidth = aisle_halfwidth, front_seat_offset = front_seat_offset) #Size for max pax count
    else
        if Vfcs_is_input == true
            if (Vfcs == 0.0 || fcs_loc == 1.0 || run_from_within == false)
                """NILS: calculate floor angle that yields maximum cabin width"""
                θ = find_floor_angles(false, Rfuse, dRfuse, h_seat = h_seat) #Find the floor angle
                pim2θ = pi - 2 * (-θ)  # NILS (see above docstring on angle sign convention)
                A_cargo = Rfuse^2 / 2 * (pim2θ - sin(pim2θ))  # NILS
                run_from_within = true  # NILS: re-run with updated θ (other conditions apply further down)
                # NILS: Afcs gets computed when run_from_within = true
            elseif (Vfcs > 0.0 && fcs_loc == 0.0 && run_from_within == true)
                """NILS: calculate floor angle that yields same cargo volume as without FCS"""
                # Afcs = Vfcs / lcyl  # NILS: required FCS cross-sectional area is volume divided by cylinder length; inaccurate for two fuselage FCS compartments separated by center wing box
                ###
                lfcs = (lcyl - wing.layout.root_chord * wing.center.cross_section.width_to_chord)  # see wingbox_coordinate_compiler.jl, and (164) and Figure 7 in tasopt.pdf
                Afcs = Vfcs / lfcs  # NILS: required FCS cross-sectional area is volume divided by usable FCS length
                ###
                pim2θ = pim2θ_from_area(Rfuse, A_cargo + Afcs)  # NILS
                θ = -(pi - pim2θ) / 2  # NILS (see above docstring on angle sign convention)
                run_from_within = false  # NILS: to avoid infinite recursion
                # NILS: A_cargo got computed when run_from_within = false
            end
        elseif Vfcs_is_input == false
            if fcs_loc == 1.0
                """NILS: calculate floor angle that yields maximum cabin width"""
                θ = find_floor_angles(false, Rfuse, dRfuse, h_seat = h_seat) #Find the floor angle
                pim2θ = pi - 2 * (-θ)  # NILS (see above docstring on angle sign convention)
                Afcs = pi * Rfuse^2 / 2  # NILS
                A_cargo = Rfuse^2 / 2 * (pim2θ - sin(pim2θ))  # NILS: area of circular segment based on https://en.wikipedia.org/wiki/Circular_segment
            elseif fcs_loc == 0.0
                """ NILS
                calculate available cross-sectional area for FCS and cargo given θ.
                See https://mit-lae.github.io/TASOPT.jl/dev/structures/cabin_sizing/#Theory-Double-decker-aircraft-7a29cdcf253b565d
                and src/utils/constants.jl for ULD dimensions. For the A320, I should use the "LD3-45" type (see
                https://www.anacargo.jp/en/int/specification/cm_container.html). See also MinCargoHeightConst(x, fuse).
                """
                ULD = fuse.cabin.unit_load_device
                ULDdims = UnitLoadDeviceDimensions[ULD]
                minwidth = ULDdims[2] #Base width
                c_ULD = minwidth  # chord of circular segment (https://en.wikipedia.org/wiki/Circular_segment)
                h_ULD = ULDdims[1]
                dh_ULD = Rfuse - sqrt(Rfuse^2 - c_ULD^2 / 4)  # sagitta of circular segment (https://en.wikipedia.org/wiki/Circular_segment)
                θfcs = -asin((Rfuse - dh_ULD - h_ULD) / Rfuse)  # angle between FCS floor and horizontal (negative downwards)
                Afcs = Rfuse^2 / 2 * ((pi + 2 * θ) - sin(pi + 2 * θ)) - (Rfuse^2 / 2 * ((pi + 2 * θfcs) - sin(pi + 2 * θfcs)))  # see 4th-Year-Scribbles on GoodNotes from 02.11.2025 (area of circular segment minus)
                Acargo = Rfuse^2 / 2 * ((pi + 2 * θfcs) - sin(pi + 2 * θfcs))  # see 4th-Year-Scribbles on GoodNotes from 02.11.2025 (area of circular segment minus)
                # println("c_ULD, h_ULD, dh_ULD, θfcs, Afcs, Acargo = $c_ULD, $h_ULD, $dh_ULD, $θfcs, $Afcs, $Acargo")  # test
            end
        end
        paxsize = fuse.cabin.exit_limit #maximum number of passengers
        w = find_cabin_width(Rfuse, wfb, nfweb, θ, h_seat) #Cabin width
        lcyl, _, seats_per_row = place_cabin_seats(paxsize, w, seat_pitch = seat_pitch, seat_width = seat_width, 
            aisle_halfwidth = aisle_halfwidth, front_seat_offset = front_seat_offset) #Cabin length
        
        # NILS: calculate underfloor FCS volume and cargo volume from cross-sectional areas and cylinder length
        # Vcargo = Acargo * lcyl  # inaccurate for two fuselage FCS compartments separated by center wing box
        lcargo = (lcyl - wing.layout.root_chord * wing.center.cross_section.width_to_chord)  # see wingbox_coordinate_compiler.jl, and (164) and Figure 7 in tasopt.pdf
        Vcargo = Acargo * lcargo
        if Vfcs_is_input == true
            nothing
        elseif Vfcs_is_input == false
            if fcs_loc == 1.0
                Vfcs = lfcs * pi * Rfuse^2  # volume of disk with radius Rfuse and length lfcs
            elseif fcs_loc == 0.0
                # Vfcs = Afcs * lcyl  # inaccurate for two fuselage FCS compartments separated by center wing box
                """NILS: calculate combined volume of two fuselage FC compartments, where
                the forward comparment starts at fuse.layout.x_start_cylinder and ends at
                wing.layout.box_x - wing.layout.root_chord * wing.center.cross_section.width_to_chord / 2
                (wing box location minus half wing box root chord - see Figure 7 in tasopt.pdf),
                and the aft compartment ends at fuse.layout.x_end_cylinder and starts at
                wing.layout.box_x + wing.layout.root_chord * wing.center.cross_section.width_to_chord / 2
                (wing box location minus half wing box root chord - see Figure 7 in tasopt.pdf).
                """
                lfcs = (lcyl - wing.layout.root_chord * wing.center.cross_section.width_to_chord)  # see wingbox_coordinate_compiler.jl, and (164) and Figure 7 in tasopt.pdf
                Vfcs = Afcs * lfcs
                lfcs_front = (wing.layout.box_x - wing.layout.root_chord * wing.center.cross_section.width_to_chord / 2) - fuse.layout.x_start_cylinder
                Vfcs_front = Afcs * lfcs_front
                lfcs_aft = fuse.layout.x_end_cylinder - (wing.layout.box_x + wing.layout.root_chord * wing.center.cross_section.width_to_chord / 2)
                Vfcs_aft = Afcs * lfcs_aft
                # Vfcs = Vfcs_front + Vfcs_aft
                println("Vfcs_front, Vfcs_aft, Vfcs_front + Vfcs_aft, Vfcs = $Vfcs_front, $Vfcs_aft,  $(Vfcs_front + Vfcs_aft), $Vfcs")
                ###
            end
        end
    end

    # NILS: log available FCS volume
    # parg[igVfcsfus] = Vfcs
    ac.nils.V_fcs_fuselage = Vfcs

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

    # NILS: re-run function to re-size cabin so as to maintain A_cargo
    if Vfcs_is_input == true
        if (Vfcs > 0.0 && fcs_loc == 0.0 && run_from_within == true)
            seats_per_row = update_fuse_for_pax!(
                ac; run_from_within=run_from_within, lcyl=lcyl, A_cargo=A_cargo,  # NILS: note the semicolon separating positional from keyword arguments
            )
        end
    elseif Vfcs_is_input == false
        nothing
    end

    #Update positions and fuselage length
    """ NILS
    It has been found that increasing lcycl itself,
    as opposed to adding `lfcs` to x_end_cylinder,
    resulted in a consistent fuselage geometry.
    """
    lcyl = lcyl + fuse.cabin.rear_seat_offset + lfcs #Make cabin longer to leave room in the back (lfcs added by NILS)

    #Update positions and fuselage length
    fuse.layout.x_end_cylinder = fuse.layout.x_start_cylinder + lcyl

    # NILS: log longitudinal location of FCS centroid
    """ NILS
    The following is correct because, despite the fact that fuse.layout.x_end_cylinder gets
    updated in update_fuse! (see below), the cabin length `lcyl` and `fuse.layout.x_start_cylinder`
    remain constant, and both fuselage FCS integration options are located before the tanks.

    #Update positions and fuselage length
    fuse.layout.x_end_cylinder = fuse.layout.x_start_cylinder + nftanks * (lftank + lftoffset) + lcyl
    """
    if fcs_loc == 1.0
        fuse.layout.x_centroid_fcs = fuse.layout.x_end_cylinder - lfcs / 2
    elseif fcs_loc == 0.0
        # fuse.layout.x_centroid_fcs = fuse.layout.x_start_cylinder + lcyl / 2  # inaccurate for two fuselage FCS compartments separated by center wing box
        """NILS: calculate combined centroid of two fuselage FC compartments, where
        the forward comparment starts at fuse.layout.x_start_cylinder and ends at
        wing.layout.box_x - wing.layout.root_chord * wing.center.cross_section.width_to_chord / 2
        (wing box location minus half wing box root chord - see Figure 7 in tasopt.pdf),
        and the aft compartment ends at fuse.layout.x_end_cylinder and starts at
        wing.layout.box_x + wing.layout.root_chord * wing.center.cross_section.width_to_chord / 2
        (wing box location minus half wing box root chord - see Figure 7 in tasopt.pdf).
        """
        x_centroid_fcs_front = fuse.layout.x_start_cylinder + ((wing.layout.box_x - wing.layout.root_chord * wing.center.cross_section.width_to_chord / 2) - fuse.layout.x_start_cylinder) / 2
        x_centroid_fcs_aft = fuse.layout.x_end_cylinder - (fuse.layout.x_end_cylinder - (wing.layout.box_x + wing.layout.root_chord * wing.center.cross_section.width_to_chord / 2)) / 2
        fuse.layout.x_centroid_fcs = (x_centroid_fcs_front * Vfcs_front + x_centroid_fcs_aft * Vfcs_aft) / Vfcs
    end

    #Update wingbox position
    wing.layout.box_x = fuse.layout.x_start_cylinder + wbox_cabin_frac * lcyl
       
    #Update other lengths
    fuse.layout.x_pressure_shell_aft = fuse.layout.x_end_cylinder + dxcyl2shellaft

    fuse.layout.x_cone_end = fuse.layout.x_pressure_shell_aft + dxshell2conend
    fuse.APU.r = [fuse.layout.x_pressure_shell_aft + dxshell2apu, 0.0,0.0]
    fuse.layout.x_end = fuse.APU.x + dxapu2end
    fuse.HPE_sys.r = [fuse.layout.x_cone_end * 0.52484, 0.0, 0.0] #TODO: address this
    
    htail.layout.box_x = fuse.layout.x_cone_end - dxhbox2conend
    vtail.layout.box_x = fuse.layout.x_cone_end - dxvbox2conend
    
    parg[igxeng    ] =  wing.layout.box_x - dxeng2wbox #Move engine

    fuse.layout.l_cabin_cylinder = lcyl #Store new cabin length

    EvaluateCabinProps!(fuse) #Update cabin parameters

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