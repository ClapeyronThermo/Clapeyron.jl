abstract type NaNLiquidModel <: LiquidVolumeModel end
@newmodelsingleton NaNLiquid NaNLiquidModel

"""
    NaNLiquid <: IdealModel

    NaNLiquid(components; 
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters
None
## Description
Placeholder liquid model. Returns NaN for any input.

## Model Construction Examples
```
# Because this model does not have parameters, all those constructors are equivalent:
liquid = NaNLiquid()
liquid = NaNLiquid("water")
liquid = NaNLiquid(["water","carbon dioxide"])
```
"""
NaNLiquid

export NaNLiquid

function volume_impl(model::NaNLiquidModel,p,T,z=SA[1.0],phase=:unknown,threaded=false,vol0 = 0.0)
    _0 = zero(first(z))
    return _0/_0
end

check_arraysize(::NaNLiquid,x::AbstractVector) = nothing
check_arraysize(::NaNLiquid,x::AbstractMatrix) = nothing
