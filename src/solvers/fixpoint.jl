abstract type AbstractFixPoint end

"""
  Solvers.fixpoint(f,x0::Real,method=SimpleFixPoint())
does a fixpoint iteration convergence on a series of real numbers.
f is a function that should 
the following strategies:
 SimpleFixPoint(): simple fixpoint, evaluates until criteria are met
 AitkenFixPoint(): uses Aitken's delta-squared process to accelerate convergence of the series
"""

struct SimpleFixPoint <: AbstractFixPoint end
struct AitkenFixPoint <: AbstractFixPoint end

function fixpoint(f,x0,
    method::AbstractFixPoint = SimpleFixPoint();
    atol=zero(float(real(first(x0)))),
    rtol=8eps(one(float(real(first(x0))))), 
    max_iters=100)
    x0,atol,retol = promote(x0,atol,rtol)
return fixpoint(f,x0,method,atol,rtol,max_iters)
end

function fixpoint(f,
        x0::T,
        ::SimpleFixPoint,
        atol::T = zero(T),
        rtol::T =8*eps(T),
        max_iters=100) where {T}

    
    nan = (0*x0)/(0*x0)
    
    xi = f(x0)
    itercount = 1 #one iteration is already done
    #if fixed point is already the initial point
    xi == x0 && return xi
    isnan(xi) && return nan

    xold = x0
    while itercount < max_iters
    xi = f(xi)
    xi == xold && return xi
    isnan(xi) && return nan
    Δx = xi-xold
    #@show Δx,xi
    if abs(Δx) < max(atol,abs(xi)*rtol)
        return xi
    end
    #println("iter = $itercount")
    itercount +=1
    xold = xi    
end
    return nan
end

function fixpoint(f,
    x0::T,
    ::AitkenFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100) where {T}


nan = (0*x0)/(0*x0)
itercount = 0
xi = x0
xii = f(xi)
xi == xii && return xii
isnan(xii) && return nan
itercount += 1

xiii = f(xii)
xii == xiii && return xiii
isnan(xiii) && return nan
@show xi, xii, xiii
itercount += 1


while itercount < max_iters
    xii = f(xi)
    xi == xii && return xii
    isnan(xii) && return nan
    itercount += 1
    
    xiii = f(xii)
    xii == xiii && return xiii
    isnan(xiii) && return nan
    @show xi, xii, xiii
    itercount += 1

    denom = (xiii - xii) - (xii - xi)
    @show denom
    xnew = xi - ((xii - xi)^2 / (xiii - 2 * xii + xi))
    #xnew = xiii - abs2(xiii-xii)/denom
    xnew == xiii && return xnew
    isnan(xnew) && return nan

    Δx = xnew-xiii
    @show Δx
    if abs(Δx) < max(atol,abs(xnew)*rtol)
        return xiii
    end
    @show xnew
    itercount +=1
    xi = xnew
end
return nan
end
#vc = 101325.9999968868
#vn   101325.99999857985