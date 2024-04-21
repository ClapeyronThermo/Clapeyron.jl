

abstract type AbstractFixPoint end

"""
  Solvers.fixpoint(f,x0::Real,method=SSFixPoint())
does a fixpoint iteration convergence on a series of real numbers.
f is a function that should 
the following strategies:
 - `SSFixPoint(dampingfactor = 1.0)`: performs succesive substitutions until convergence is met. 
α = `dampingfactor` is determines a buffer for each iteration, defined as `x1 = α*f(x1) + (1-α)*x1`

- `Aitken()`: uses Aitken's delta-squared process to accelerate convergence of the series. recommended for harmonic iterates.

"""
function fixpoint end

struct SSFixPoint{T<:Real} <: AbstractFixPoint 
    dampingfactor::T
    lognorm::Bool
    normorder::Float64
end

SSFixPoint(;dampingfactor=1.0,lognorm = false,normorder = Inf) = SSFixPoint(dampingfactor,lognorm,Float64(normorder))
SSFixPoint(dampingfactor) = SSFixPoint(dampingfactor,false,Inf)
function promote_method(method::SSFixPoint,T)
    return SSFixPoint(T(method.dampingfactor),method.lognorm,method.normorder)
end

struct AitkenFixPoint <: AbstractFixPoint end

function promote_method(method::AitkenFixPoint,T)
    return method
end

function convergence(xold,xi,atol,rtol,lognorm = false,normorder = Inf)
    not_finite = false
    for xii in xi
        if !isfinite(xii)
            not_finite = true
            break
        end
    end
    not_finite && return (true,false) #terminate, with nan
    xi == xold && return (true,true) #terminate, with current number
    if xi isa Number
        if lognorm
            Δx = abs(xi-xold)
        else
            Δx = abs(xi/xold - 1)
        end
    else
        if lognorm
            Δx = norm(((xi[i]/xold[i] - 1) for i in eachindex(xold,xi)),normorder)
        else
            Δx = norm((xi[i] - xold[i] for i in eachindex(xold,xi)),normorder)
        end
    end
    if lognorm
        #ignore rtol in lognorm
        normxi = zero(eltype(xi))/one(eltype(xi))
    else
        normxi = norm(xi,normorder)
    end

    if abs(Δx) < max(atol,normxi*rtol)
        return (true,true) #terminate, with current number
    end
    return (false,false) #keep iterating
end

function fixpoint(f,x0,
    method::AbstractFixPoint = SSFixPoint();
    atol=zero(eltype(x0)),
    rtol=8eps(one(eltype(x0))), 
    max_iters=100,
    return_last = false)
    _,atol,rtol = promote(one(eltype(x0)),atol,rtol)
    method = promote_method(method,eltype(x0))
    return _fixpoint(f,x0,method,atol,rtol,max_iters,return_last)
end

function _fixpoint(f::F,
    x0::T,
    method::SSFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100,
    return_last = false) where {F,T<:Real}
    
    nan = (0*atol)/(0*atol)
    xi = f(x0)
    α = method.dampingfactor
    ℕ = method.normorder
    lognorm = method.lognorm
    converged,finite = convergence(x0,xi,atol,rtol,lognorm,ℕ)
    converged && return ifelse(finite,xi,nan)
    itercount = 1
    xold = x0
    while itercount < max_iters
        xi = α*f(xi) + (1-α)*xi  
        converged,finite = convergence(xold,xi,atol,rtol,lognorm,ℕ)
        converged && return ifelse(finite,xi,nan)    
        itercount +=1
        xold = xi
    end
    return ifelse(return_last,xi,nan)
end

function _fixpoint(f::F,
    x00::T,
    method::AitkenFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100,
    return_last = false) where {F,T<:Real}

    nan = (0*atol)/(0*atol)
    itercount = 2
    x3 = x00
    x1 = f(x3)
    converged,finite = convergence(x3,x1,atol,rtol)
    converged && return ifelse(finite,x1,nan)  
    x2 = f(x1)
    converged,finite = convergence(x1,x2,atol,rtol)
    converged && return ifelse(finite,x2,nan)

    while itercount < max_iters
        itercount += 1
        λ2 = (x2 - x1)/(x1-x3)
        dx = -(λ2/(1-λ2))*(x2-x1)
        x3 = x2 + dx
        converged,finite = convergence(x2,x3,atol,rtol)
        converged && return ifelse(finite,x3,nan)
        
        itercount += 1
        x1 = f(x3)
        converged,finite = convergence(x3,x1,atol,rtol)
        converged && return ifelse(finite,x1,nan)
        
        itercount += 1
        x2 = f(x1)
        converged,finite = convergence(x1,x2,atol,rtol)
        converged && return ifelse(finite,x2,nan)
    end
    return ifelse(return_last,x2,nan)
end

function _fixpoint(f!::F,
    x0::X where {X <:AbstractVector{T}},
    method::SSFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100,
    return_last = false) where {F,T<:Real}
    
    nan = (0*atol)/(0*atol)
    xi = copy(x0)
    xi = f!(xi,x0)
    α = method.dampingfactor
    ℕ = method.normorder
    lognorm = method.lognorm
    converged,finite = convergence(x0,xi,atol,rtol,lognorm,ℕ)
    if converged
        if finite
            return xi
        else
            xi .= nan
            return xi 
        end
    end
    itercount = 1
    xold = copy(x0)
    while itercount < max_iters
        xi = f!(xi,xold)
        xi .*= α
        xi .+= (1 .- α) .* xold
        converged,finite = convergence(xold,xi,atol,rtol,lognorm,ℕ)
        if converged
            if finite
                return xi
            else
                xi .= nan
                return xi 
            end
        end
        itercount +=1
        xold .= xi
    end
    !return_last && (xi .= nan)
    return xi
end
