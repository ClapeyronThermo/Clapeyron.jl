struct BasicIdealParam <: EoSParam
end

abstract type BasicIdealModel <: IdealModel end
struct BasicIdeal <: BasicIdealModel
    params::BasicIdealParam
end

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

export BasicIdeal
function BasicIdeal(components::Vector; userlocations::Array{String,1}=String[], verbose=false)
    return BasicIdeal(BasicIdealParam())
end


is_splittable(::BasicIdeal) = false
recombine_impl!(model::BasicIdeal) = model
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