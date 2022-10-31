abstract type BasicIdealModel <: IdealModel end

BasicIdeal_SETUP = ModelOptions(
        :BasicIdeal;
        supertype=BasicIdealModel,
        has_components=false,
        has_params=false,
    )

createmodel(BasicIdeal_SETUP; verbose=true)
export BasicIdeal

"""
    BasicIdeal <: IdealModel
    BasicIdeal(components::Array{String,1}; 
    userlocations::Array{String,1}=String[], 
    verbose=false)

## Input parameters

None

## Description

Default Ideal Model. Constant specific heat capacity equal to `5R/2`. it's Helmholtz energy is equal to:
```
    a₀ = A₀/nRT =  ∑(xᵢlog(nxᵢ/V)) - 1 - 1.5log(T)
```
"""
BasicIdeal

is_splittable(::BasicIdeal) = false

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
