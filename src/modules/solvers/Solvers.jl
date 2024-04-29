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

function st_solve(B,∇f,::Val{1})
    ∇f̄ = SVector(f[1])
    B̄ = SMatrix{1,1}((B[1,1],))
    return B̄\∇f
end

function st_solve(B,∇f,::Val{2})
    ∇f̄ = SVector(∇f[1],∇f[2])
    B̄ = SMatrix{2,2}((B[1,1],B[2,1],B[1,2],B[2,2]))
    return B̄\∇f
end

function st_solve(B,∇f,::Val{3})
    ∇f̄ = SVector(∇f[1],∇f[2],∇f[3])
    B̄ = SMatrix{3,3}((B[1,1],B[2,1],B[3,1],B[1,2],B[2,2],B[3,2],B[1,3],B[2,3],B[3,3]))
    return B̄\∇f
end

function st_solve(B,∇f,::Val{4})
    ∇f̄ = SVector(∇f[1],∇f[2],∇f[3],∇f[4])
    B̄ =  SMatrix{4,4}((B[1,1],B[2,1],B[3,1],B[4,1],B[1,2],B[2,2],B[3,2],B[4,2],B[1,3],B[2,3],B[3,3],B[4,3],B[1,4],B[2,4],B[3,4],B[4,4]))
    return B̄\∇f
end


function try_st_linsolve(d,B,∇f)
    n = length(∇f)
    success = 1 <= n <= 3
    if n == 1
        d .= st_solve(B,∇f,Val{1}())
    elseif n == 2
        d .= st_solve(B,∇f,Val{2}())
    elseif n == 3
        d .= st_solve(B,∇f,Val{3}())
    elseif n == 4
        d .= st_solve(B,∇f,Val{4}())
    end
    return d,success
end

function try_st_linsolve(B,∇f)
    d = similar(∇f)
    return try_st_linsolve(d,B,∇f)
end

function cholesky_linsolve(d,B,∇f)
    d,success = try_st_linsolve(d,B,∇f)
    success && return d
    cholesky!(Positive, B)
    Bchol = Cholesky(B,'L',0)
    d .=  Bchol\∇f
end

Base.summary(::NLSolvers.Newton{<:Direct, typeof(cholesky_linsolve)}) = "Newton's method with Cholesky linsolve"
CholeskyNewton() = NLSolvers.Newton(linsolve=cholesky_linsolve)

function static_linsolve(d,B,∇f)
    d,success = try_st_linsolve(d,B,∇f)
    success && return d
    d .=  B\∇f
end

function static_linsolve(B,∇f)
    if length(∇f) <= 3
    d,_ = try_st_linsolve(B,∇f)
        return d
    else
        return B\∇f
    end
end

Base.summary(::NLSolvers.Newton{<:Direct, typeof(static_linsolve)}) = "Newton's method with optional static linsolve"


const LUPivot = @static if VERSION < v"1.7beta"
    Val(true)
else
    RowMaximum()
end

function lup_linsolve(d,B,∇f)
    d,success = try_st_linsolve(d,B,∇f)
    success && return d
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

export CholeskyNewton,LUPNewton,static_linsolve

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
include("fixpoint/anderson.jl")
include("optimize.jl")
include("integral.jl")
include("chebyshev.jl")


end # module
