abstract type BasicIdealModel <: IdealModel end
@newmodelsingleton BasicIdeal BasicIdealModel

"""
    BasicIdeal <: IdealModel
    BasicIdeal(components::Array{String,1}; 
    userlocations::Array{String,1}=String[], 
    verbose = false)
## Input parameters
None
## Description
Default Ideal Model. Constant specific heat capacity equal to `5R/2`. it's Helmholtz energy is equal to:
```
    a₀ = A₀/nRT =  ∑(xᵢlog(nxᵢ/V)) - 1 - 1.5log(T)
```

## Model Construction Examples
```
# Because this model does not have parameters, all those constructors are equivalent:
idealmodel = BasicIdeal()
idealmodel = BasicIdeal("water")
idealmodel = BasicIdeal(["water","carbon dioxide"])
```
"""
BasicIdeal

export BasicIdeal

function a_ideal(model::BasicIdeal, V, T, z)
    N = ∑(z)
    #x = z/∑(z)
    res = ∑(xlogx,z) 
    res /= N 
    res -= log(V) 
    res -= 1.5*log(T)
    res -= one(res)
    # ∑(x .* log.(z/V)) - 1 original formulation, prone no NaN when passing pure Fractions
    return res
end

check_arraysize(::BasicIdealModel,x::AbstractVector) = nothing
check_arraysize(::BasicIdealModel,x::AbstractMatrix) = nothing
