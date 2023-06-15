"""
    calculate_missing_mixing(params::EmpiricMultiFluidParam,mixing::MixingRule)

Calculates missing values, using the parameters stored in `params`. modifies `mixing` implace. this function is called at `EmpiricMultiFluid` model creation time.

"""
recombine_mixing!(model::EmpiricMultiFluid,mixing) = nothing



include("Asymmetric.jl")
include("LB.jl")
include("linear.jl")
include("TL.jl")

v_scale(model::EmpiricMultiFluid,z,mixing::Nothing,∑z = sum(z)) = v_scale(model,z,LinearMixing(),∑z)
T_scale(model::EmpiricMultiFluid,z,mixing::Nothing,∑z = sum(z)) = T_scale(model,z,LinearMixing(),∑z)

