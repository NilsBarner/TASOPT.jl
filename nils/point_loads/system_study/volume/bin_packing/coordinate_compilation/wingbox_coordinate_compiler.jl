""" This scrip calculate the x-y-z coordinates of the wing box
for plotting and bin packing purposes. It is based on stickfig()
in src/IO/output_plots.jl and scales the wing planform to the
wing box planform as described in the below code."""

export compile_wingbox_coordinates

function compile_wingbox_coordinates(ac::aircraft)

    #if aircraft is not sized, cannot plot
    if !ac.is_sized[1] 
        @warn "The aircraft ($(ac.name)) must be sized before being plotted. Skipping `stick_fig`..."
        return nothing
    end

    # Unpack aircraft components
    parg, parm, para, pare, options, fuse, fuse_tank, wing, htail, vtail, engine = unpack_ac(ac,1) #imission = 1 

    # Wing
    co = wing.layout.root_chord
    cs = wing.layout.root_chord*wing.inboard.λ
    ct = wing.layout.root_chord*wing.outboard.λ

    # NILS: Go from wing chords to wing box chords
    # NOTE that, while wing.center.cross_section.width_to_chord is defined
    # perpendicular to flow, for the scaling from wing to wing box chord
    # it does not matter as both get scaled by the same factor cosd(wing.layout.sweep)
    # See (164) and Figure 7 in tasopt.pdf, and calc_sparbox_internal_area() in wing_weights.jl.
    co *= wing.center.cross_section.width_to_chord
    cs *= wing.inboard.cross_section.width_to_chord
    ct *= wing.outboard.cross_section.width_to_chord

    sweep = wing.layout.sweep
    λs = wing.inboard.λ
    λt = wing.outboard.λ

    bo = wing.layout.root_span
    bs = wing.layout.break_span
    b  = wing.layout.span

    # NILS: go from wing thickness to wing box thickness
    # I will here make the simplifying assumption of a
    # rectangular cross-section as opposed to one where the
    # "[...] height is assumed to taper off quadratically to
    # a fraction rh at the webs, [...]" (page 34 in tasopt.pdf).
    web_height_s = wing.inboard.cross_section.web_to_box_height * wing.inboard.cross_section.thickness_to_chord * co

    xax = 0.40
    xcLE = -xax
    xcTE = 1.0 - xax

    dx = wing.layout.box_x    
    etas = bs/b
    etao = bo/b
    cosL = cos(sweep*pi/180.0)
    tanL = tan(sweep*pi/180.0)

    xs = tanL*(bs-bo)/2.0
    xt = tanL*(b -bo)/2.0

    xw = zeros(12)
    yw = zeros(12)
    zw = zeros(12)

    #X locations of wing vertices
    xw[1] =      co*xcLE + dx
    xw[2] = xs + cs*xcLE + dx
    xw[3] = xt + ct*xcLE + dx
    xw[4] = xt + ct*xcTE + dx
    xw[5] = xs + cs*xcTE + dx
    xw[6] =      co*xcTE + dx

    xw[7] = xw[1]
    xw[8] = xw[2]
    xw[9] = xw[3]
    xw[10] = xw[4]
    xw[11] = xw[5]
    xw[12] = xw[6]
    
    #Y locations of wing vertices
    yw[1] = bo/2.0
    yw[2] = bs/2.0
    yw[3] = b /2.0
    yw[4] = b /2.0
    yw[5] = bs/2.0
    yw[6] = bo/2.0

    yw[7] = yw[1]
    yw[8] = yw[2]
    yw[9] = yw[3]
    yw[10] = yw[4]
    yw[11] = yw[5]
    yw[12] = yw[6]

    #Z locations of wing vertices (added by NILS)
    zw[1] = 0
    zw[2] = 0
    zw[3] = 0
    zw[4] = 0
    zw[5] = 0
    zw[6] = 0

    zw[7] = web_height_s
    zw[8] = web_height_s * cs / co
    zw[9] = web_height_s * ct / co
    zw[10] = web_height_s * ct / co
    zw[11] = web_height_s * cs / co
    zw[12] = web_height_s

    return xw, yw, zw

end