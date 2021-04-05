

abstract type AbstractFixPoint end

"""
  Solvers.fixpoint(f,x0::Real,method=SimpleFixPoint())
does a fixpoint iteration convergence on a series of real numbers.
f is a function that should 
the following strategies:
 SimpleFixPoint(): simple fixpoint, evaluates until criteria are met
 AitkenFixPoint(): uses Aitken's delta-squared process to accelerate convergence of the series
"""

struct SimpleFixPoint{T<:Real} <: AbstractFixPoint 
    dampingfactor::T
end

SimpleFixPoint(;dampingfactor=1.0) = SimpleFixPoint(dampingfactor)


struct AitkenFixPoint <: AbstractFixPoint end

function fixpoint(f,x0,
    method::AbstractFixPoint = SimpleFixPoint();
    atol=zero(float(real(nested_eltype(x0)))),
    rtol=8eps(one(float(real(nested_eltype(x0))))), 
    max_iters=100,
    norm = norm)
    x0,atol,retol = promote(x0,atol,rtol)
    return _fixpoint(f,x0,method,atol,rtol,max_iters,norm)
end

function _fixpoint(f,
        x0::T,
        method::SimpleFixPoint = SimpleFixPoint(),
        atol::T = zero(T),
        rtol::T =8*eps(T),
        max_iters=100,
        norm=norm) where {T}

    
    nan = (0*atol)/(0*atol)
    
    xi = f(x0)
    itercount = 1 #one iteration is already done
    #if fixed point is already the initial point
    xi == x0 && return xi
    isnan(xi) && return nan

    xold = x0
    while itercount < max_iters
    #("iter = $itercount")
    xi = f(xi)
    xi == xold && return xi
    isnan(xi) && return nan
    Δx = abs(xi-xold)
    #@show Δx/xi,xi
    if abs(Δx) < max(atol,abs(xi)*rtol)
        return xi
    end
    
    itercount +=1
    xold = xi    
end
    return nan
end
