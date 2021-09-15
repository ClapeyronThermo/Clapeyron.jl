

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
end

SSFixPoint(;dampingfactor=1.0) = SSFixPoint(dampingfactor)

function promote_method(method::SSFixPoint,T)
    return SSFixPoint(T(method.dampingfactor))
end

struct AitkenFixPoint <: AbstractFixPoint end

function promote_method(method::AitkenFixPoint,T)
    return method
end

@inline function convergence(xold,xi,atol,rtol)
    !isfinite(xi) && return (true,false) #terminate, with nan
    xi == xold && return (true,true) #terminate, with current number
    Δx = abs(xi-xold)
    if abs(Δx) < max(atol,abs(xi)*rtol)
        return (true,true) #terminate, with current number
    end
    return (false,false) #keep iterating
end

function fixpoint(f,x0::Real,
    method::AbstractFixPoint = SSFixPoint();
    atol=zero(nested_eltype(x0)),
    rtol=8eps(one(nested_eltype(x0))), 
    max_iters=100)
    x0,atol,rtol = promote(x0,atol,rtol)
    method = promote_method(method,typeof(x0))
    return _fixpoint(f,x0,method,atol,rtol,max_iters)
end

function fixpoint(f,x0::AbstractVector{<:Real},
    method::AbstractFixPoint = AndersonFixPoint();
    atol=zero(nested_eltype(x0)),
    rtol=8eps(one(nested_eltype(x0))), 
    max_iters=100)
    x0,atol,rtol = promote(x0,atol,rtol)
    method = promote_method(method,nested_eltype(x0))
    return _fixpoint(f,x0,method,atol,rtol,max_iters)
end

function _fixpoint(f::F,
    x0::T,
    method::SSFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100) where {F,T}
    nan = (0*atol)/(0*atol)
    xi = f(x0)
    converged,finite = convergence(x0,xi,atol,rtol)
    converged && return ifelse(finite,xi,nan)
    itercount = 1
    xold = x0
    while itercount < max_iters
        xi = f(xi)  
        converged,finite = convergence(xold,xi,atol,rtol)
        converged && return ifelse(finite,xi,nan)    
        itercount +=1
        xold = xi
    end
    return nan
end

function _fixpoint(f::F,
    x00::T,
    method::AitkenFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100) where {F,T}

    α = method.dampingfactor
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
        x3 = x2 + α*dx
    
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
    return nan
end

Base.@kwdef struct AndersonFixPoint{T} <: AbstractFixPoint
    theta::T = 0.01
    tau::T = 0.001
    D::T = 1e6
    eps::T = 1e-6
    m::Int = 5
    dampingfactor::T = 0.1
end

function _fixpoint(f,
    x0::AbstractVector{<:T},
    method::AndersonFixPoint,
    atol::T = zero(T),
    rtol::T =8*eps(T),
    max_iters=100) where {T <: Real}
    
    nan = (0*atol)/(0*atol)
    θ =  method.theta
    τ = method.tau
    α = method.dampingfactor
    ϵ = method.eps
    D = method.D
    m = method.m
    xᵢ₋₁ = x0
    mᵢ = 0
    naa = 0
    xᵢ = f(xᵢ₋₁)
    Hk = one(xi)
    gi = abs(xᵢ₋₁ - xᵢ)
    U = gi
    yᵢ = gi
    gold = zero(xi)
end

