module EoSFunctions
using ..Solvers
using ForwardDiff
import LogExpFunctions


"""
    bmcs_hs(ζ0,ζ1,ζ2,ζ3)

Boublík-Mansoori-Leeland-Carnahan-Starling hard sphere term:
```
1/ζ₀ * (3ζ₁*ζ₂/(1-ζ₃) + ζ₂^3/(ζ₃*(1-ζ₃)^2) + (ζ₂^3/ζ₃^2-ζ₀)*log(1-ζ₃))
```
"""
function bmcs_hs(ζ0,ζ1,ζ2,ζ3)
    ζ3m1 = (1-ζ3)
    ζ3m1² = ζ3m1*ζ3m1
    return 1/ζ0 * (3ζ1*ζ2/ζ3m1 + ζ2^3/(ζ3*ζ3m1²) + (ζ2^3/ζ3^2-ζ0)*Solvers.log1p(-ζ3))
end

"""
    xlogx(x::Real,k = one(x))
Return `x * log(k*x)` for `x ≥ 0`, handling ``x = 0`` by taking the downward limit.

copied from LogExpFunctions.jl
"""
function xlogx(x::Real,k = one(x))
    _0 = zero(x)
    iszero(x) && return _0
    ifelse(x > _0,x*Base.log(max(_0,k*x)),_0/_0)
end

logabssinh(x) = LogExpFunctions.logabssinh(x)
logcosh(x) = LogExpFunctions.logcosh(x)

#=function logabssinh(d::ForwardDiff.Dual{T}) where {T}
    x = logabssinh(ForwardDiff.value(d))
    dx = coth(x)
    return ForwardDiff.Dual{T}(x, dx * ForwardDiff.partials(d))
end

function logcosh(d::ForwardDiff.Dual{T}) where {T}
    x = logcosh(ForwardDiff.value(d))
    dx = tanh(x)
    return ForwardDiff.Dual{T}(x, dx * ForwardDiff.partials(d))
end
=#
export bmcs_hs,xlogx

end #module
