"""
    $TYPEDEF

Field of an [`aircraft`](@ref) containing custom parameters introduced by Nils.

$TYPEDFIELDS
"""
@kwdef mutable struct Nils
    theta_floor::Float64
    engine_point_load::Float64
    TR_scale::Float64
    tdivc_scale::Float64
    V_fcs_nacelle::Float64
    sigma_fcs_nacelle::Float64
    sigma_fcs::Float64
    wing_frac::Float64
    nacelle_frac::Float64
    fcs_loc::Any
    span_loc::Any
    V_fcs_fuselage::Float64
    l_fcs_fuselage::Float64
    fcs_fuselage_location::Any
    V_fcs_wing::Float64
end

function Base.summary(nils::Nils)
  println("TO BE WRITTEN - SEE options.jl!")
end
function Base.show(io::IO, nils::Nils)
    println("TO BE WRITTEN - SEE options.jl!")
end