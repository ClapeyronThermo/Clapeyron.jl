module Solvers

using LinearAlgebra
using NLSolvers,Roots
using PositiveFactorizations
using DiffResults, ForwardDiff
using StaticArrays

"""
    dnorm(x,y,p)

Equivalent to `norm((xi-yi for (xi, yi) in zip(x, y)), p)`
"""
function dnorm(x,y,p = 2)
    return norm((xi-yi for (xi, yi) in zip(x, y)), p)
end


function cholesky_linsolve(d,B,∇f)
    cholesky!(Positive, B)
    Bchol = Cholesky(B,'L',0)
    d .=  Bchol\∇f
end


Base.summary(::NLSolvers.Newton{<:Direct, typeof(cholesky_linsolve)}) = "Newton's method with Cholesky linsolve"
CholeskyNewton() = NLSolvers.Newton(linsolve=cholesky_linsolve)

const LUPivot = @static if VERSION < v"1.7beta"
    Val(true)
else
    RowMaximum()
end

function lup_linsolve(d,B,∇f)
    F = lu!(B,LUPivot)
    ldiv!(d,F,∇f)
    d
end

function lup_linsolve(B,∇f)
    F = lu!(B,LUPivot)
    F\∇f
end

Base.summary(::NLSolvers.Newton{<:Direct, typeof(lup_linsolve)}) = "Newton's method with pivoted LU linsolve"
LUPNewton() = NLSolvers.Newton(linsolve=lup_linsolve)

export CholeskyNewton,LUPNewton

"""
    x_sol(res::NLSolvers.ConvergenceInfo)
    
Returns the scalar or vector x that solves the system of equations or is the minimizer of an optimization procedure.
"""
x_sol(res) = NLSolvers.solution(res)
function x_sol(res::NLSolvers.ConvergenceInfo{NLSolvers.BrentMin{Float64}})
    return res.info.x
end
include("poly.jl")
include("nanmath.jl")
include("ad.jl")
include("nlsolve.jl")
include("fixpoint/fixpoint.jl")
include("fixpoint/ADNewton.jl")
include("optimize.jl")
include("integral.jl")
include("chebyshev.jl")


end # module
