"""
    recombine_mixing!(model::MultiFluid,mixing::MixingRule)

Calculates missing values, using the parameters stored in `model`. modifies `mixing` implace. this function is called at `MultiFluid` model creation time.

"""
recombine_mixing!(model::MultiFluid,mixing) = nothing


include("Asymmetric.jl")
include("LB.jl")
include("linear.jl")
include("TL.jl")

v_scale(model::MultiFluid,z,mixing::Nothing,∑z = sum(z)) = v_scale(model,z,LinearMixing(),∑z)
T_scale(model::MultiFluid,z,mixing::Nothing,∑z = sum(z)) = T_scale(model,z,LinearMixing(),∑z)

