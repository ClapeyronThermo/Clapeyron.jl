"""
    __gas_model(model::EoSModel)

internal function.
provides the model used to calculate gas properties.
normally, this is the identity, but `CompositeModel` has a gas model by itself.
"""
__gas_model(model::EoSModel) = model

include("bubble_point.jl")
include("dew_point.jl")
include("LLE_point.jl")