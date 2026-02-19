abstract type ZeroLiquidModel <: LiquidVolumeModel end
@newmodelsingleton ZeroLiquid ZeroLiquidModel

"""
    ZeroLiquid <: IdealModel

    ZeroLiquid(components; 
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters
None
## Description
Placeholder liquid model. Returns 0 for any input.

## Model Construction Examples
```
# Because this model does not have parameters, all those constructors are equivalent:
liquid = ZeroLiquid()
liquid = ZeroLiquid("water")
liquid = ZeroLiquid(["water","carbon dioxide"])
```
"""
ZeroLiquid

export ZeroLiquid

@deprecate NaNLiquid(args...;kwargs...) ZeroLiquid(args...;kwargs...)

function volume_impl(model::ZeroLiquidModel,p,T,z,phase,threaded,vol0)
    _0 = zero(first(z))
    return _0
end

check_arraysize(::ZeroLiquidModel,x::AbstractVector) = nothing
check_arraysize(::ZeroLiquidModel,x::AbstractMatrix) = nothing
