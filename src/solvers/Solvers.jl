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
include("integral21.jl")

end # module
