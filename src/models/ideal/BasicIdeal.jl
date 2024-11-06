abstract type BasicIdealModel <: IdealModel end
@newmodelsingleton BasicIdeal BasicIdealModel

"""
    BasicIdeal <: IdealModel

    BasicIdeal(components;
    userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Input parameters
None
## Description
Default Ideal Model. Constant specific heat capacity equal to `5R/2`. it's Helmholtz energy is equal to:
```
    a₀ = A₀/nRT = ∑(xᵢlog(nxᵢ/V)) - 1 - 1.5log(T)
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
    N = sum(z)
    res = sum(Base.Fix2(xlogx,1/V),z)
    res /= N
    res -= 1.5*log(T)
    res -= one(res)
    # ∑(x .* log.(z/V)) - 1 original formulation, prone no NaN when passing pure Fractions
    return res
end

function a_ideal_T(model::BasicIdeal,T,z)
    return -1.5*log(T) - 1
end

check_arraysize(::BasicIdealModel,x::AbstractVector) = nothing
check_arraysize(::BasicIdealModel,x::AbstractMatrix) = nothing

function ∂²f∂T²(model::BasicIdeal,V,T,z)
    return -1.5*sum(z)*Rgas(model)/T
end

@newmodelsingleton ZeroIdeal BasicIdealModel

function a_ideal(model::ZeroIdeal, V, T, z)
    return zero(Base.promote_eltype(T,z))
end
#=
function ∂²f∂T²(model::CpIdeal,V,T,z)
    return -(75.34 - 8.314)*sum(z)/T
end =#
#=
function ∂²f∂T²(model::CpIdeal,V,T,z)
    return -75.34*sum(z)*Rgas(model)/T
end =#