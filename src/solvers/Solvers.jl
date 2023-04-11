module Solvers

using LinearAlgebra
using NLSolvers,Roots
using PositiveFactorizations
using DiffResults, ForwardDiff
using StaticArrays

function cholesky_linsolve(d,B,∇f)
    cholesky!(Positive, B)
    Bchol = Cholesky(B,'L',0)
    d .=  Bchol\∇f
end

function cholesky_linsolve(B,∇f)
    Bchol =cholesky(Positive, B)
    Bchol\∇f
end

Base.summary(::NLSolvers.Newton{<:Direct, typeof(cholesky_linsolve)}) = "Newton's method with Cholesky linsolve"
CholeskyNewton() = NLSolvers.Newton(linsolve=cholesky_linsolve)
export CholeskyNewton

"""
    x_sol(res::NLSolvers.ConvergenceInfo)
    
Returns the scalar or vector x that solves the system of equations or is the minimizer of an optimization procedure.
"""
x_sol(res) = NLSolvers.solution(res)

include("poly.jl")
include("ad.jl")
include("nanmath.jl")
include("nlsolve.jl")
include("fixpoint/fixpoint.jl")
include("fixpoint/ADNewton.jl")
include("optimize.jl")
include("integral.jl")


# Order 3 bracketing method
#https://github.com/JuliaMath/Roots.jl/issues/359
using Roots
import Roots: @set!, bracket, adjust_bracket, NullTracks, isissue, log_message
struct ChebyshevBracket <: Roots.AbstractBracketingMethod end
Roots.fn_argout(::ChebyshevBracket) = 3
struct ChebyshevBracketState{T,S} <: Roots.AbstractUnivariateZeroState{T,S}
    xn1::T
    xn0::T
    fxn1::S
    fxn0::S
end

# use xatol, xrtol only, but give some breathing room over the strict ones and cap number of steps
function Roots.default_tolerances(::ChebyshevBracket, ::Type{T}, ::Type{S}) where {T,S}
    xatol = eps(real(T)) * oneunit(real(T))
    xrtol = 2eps(real(T))  # unitless
    atol = 4 * eps(real(float(S))) * oneunit(real(S))
    rtol = 4 * eps(real(float(S))) * one(real(S))
    maxevals = 40
    strict = false
    (xatol, xrtol, atol, rtol, maxevals, strict)
end


function Roots.init_state(M::ChebyshevBracket, F::Roots.Callable_Function, x)
    x₀, x₁ = adjust_bracket(x)
    fx₀, fx₁ = first(F(x₀)), first(F(x₁))
    state = Roots.init_state(M, F, x₀, x₁, fx₀, fx₁)
end

function Roots.init_state(L::ChebyshevBracket, F, x₀, x₁, fx₀, fx₁)
    # keep xn1 the smallest of fxn1, fxn0
    if abs(fx₀) > abs(fx₁)
        state = ChebyshevBracketState(x₁, x₀, fx₁, fx₀)
    else
        state = ChebyshevBracketState(x₀, x₁, fx₀, fx₁)
    end
    state
end

# bracketed Chebyshev
function Roots.update_state(
    M::ChebyshevBracket,
    F,
    o::ChebyshevBracketState{T,S},
    options,
    l=NullTracks(),
) where {T,S}

    c,b,fc,fb = o.xn1, o.xn0, o.fxn1, o.fxn0
    if c < b
        a, fa = c, fc
    else
        a,b,fa,fc = b,c,fb,fc
    end

    fc::S, (Δ, Δ₂) = F(c)

    if isissue(Δ)
        log_message(l, "Issue with computing `Δ`")
        return (o, true)
    end

    # try chebyshev
    # if issue, then newton
    # if issue, then secant
    # if issue, then bisection
    x = c - Δ - 1/2 * Δ^2/Δ₂
    if !(a < x < b)
        x = c - Δ
    end

    if a < x < b
        fx,_ = F(x)
        ā,b̄,c,fā,fb̄,fc = bracket(a,b,x,fa,fb,fx)
    else
        # secant step
        x = a - fa * (b-a)/(fb-fa)
        if x - a < 0.0001 * (b-a) || b-x < 0.0001 * (b-a)
            x = a + (b-a)*0.5
        end

        fx,_ = F(x)
        ā,b̄,c,fā,fb̄,fc = bracket(a,b,x,fa,fb,fx)
    end
    c, fc = abs(fā) < abs(fb̄) ? (ā,fā) : (b̄, fb̄)

    a,b,fa,fb = ā, b̄, fā, fb̄

    @set! o.xn0 = a == c ? b : a
    @set! o.xn1 = c
    @set! o.fxn0 = a == c ? fb : fa
    @set! o.fxn1 = fc

    return o, false
end

export ChebyshevBracket

end # module
