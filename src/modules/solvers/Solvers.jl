module Solvers

using LinearAlgebra
using LinearAlgebra: BlasInt
using NLSolvers,Roots
using PositiveFactorizations
using DiffResults, ForwardDiff
using StaticArrays

struct ΔVector{T,V1,V2} <: AbstractVector{T}
    v1::V1
    v2::V2
    log::Bool
end

function ΔVector(v1::AbstractVector{T},v2::AbstractVector{T},log = false) where T
    V1,V2 = typeof(v1),typeof(v2)
    return ΔVector{T,V1,V2}(v1,v2,log)
end

Base.size(v::ΔVector) = size(v.v1)

Base.@propagate_inbounds function Base.getindex(v::ΔVector{T},i::Int) where T
    v1,v2 = v.v1[i],v.v2[i]
    if v.log
        return T(v1/v2 - 1)
    else
        return T(v1 - v2)
    end
end

"""
    dnorm(x,y,p)

Equivalent to `norm((xi-yi for (xi, yi) in zip(x, y)), p)`
"""
function dnorm(x,y,p = 2)
    Δ = ΔVector(x,y)
    return norm(Δ, p)
end

function dnorm(x::Number,y::Number,p = 2)
    return norm(x - y, p)
end

function st_solve(B,∇f,::Val{1})
    ∇f̄ = SVector(∇f[1])
    B̄ = SMatrix{1,1}((B[1,1],))
    return B̄\∇f
end

function st_solve(B,∇f,::Val{2})
    ∇f̄ = SVector(∇f[1],∇f[2])
    B̄ = SMatrix{2,2}((B[1,1],B[2,1],B[1,2],B[2,2]))
    return inv(B̄)*∇f
end

function st_solve(B,∇f,::Val{3})
    ∇f̄ = SVector(∇f[1],∇f[2],∇f[3])
    B̄ = SMatrix{3,3}((B[1,1],B[2,1],B[3,1],B[1,2],B[2,2],B[3,2],B[1,3],B[2,3],B[3,3]))
    return inv(B̄)*∇f
end

function st_solve(B,∇f,::Val{4})
    ∇f̄ = SVector(∇f[1],∇f[2],∇f[3],∇f[4])
    B̄ =  SMatrix{4,4}((B[1,1],B[2,1],B[3,1],B[4,1],B[1,2],B[2,2],B[3,2],B[4,2],B[1,3],B[2,3],B[3,3],B[4,3],B[1,4],B[2,4],B[3,4],B[4,4]))
    return inv(B̄)*∇f
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
    lup_linsolve(d,B,∇f)
end

function static_linsolve(B,∇f)
    if length(∇f) <= 3
        d,_ = try_st_linsolve(B,∇f)
        return d
    else
        return lup_linsolve(B,∇f)
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


"""
    x_sol(res::NLSolvers.ConvergenceInfo)
    
Returns the scalar or vector x that solves the system of equations or is the minimizer of an optimization procedure.
"""
x_sol(res) = NLSolvers.solution(res)
function x_sol(res::NLSolvers.ConvergenceInfo{NLSolvers.BrentMin{Float64}})
    return res.info.x
end

struct RestrictedLineSearch{F,LS} <: NLSolvers.LineSearcher
    f::F #function that restricts the line search
    ls::LS #actual line search
end

#=
unsafe LU
it does not check, it does not allow to select a pivot strategy
it allows to pass a buffer ipiv.
=#
function unsafe_LU!(A::AbstractMatrix{T}, ipiv = Vector{BlasInt}(undef, min(size(A)...))) where {T}
    #check && LAPACK.chkfinite(A)
    # Extract values
    m, n = size(A)
    minmn = min(m,n)

    # Initialize variables
    info = 0
    #ipiv = Vector{BlasInt}(undef, minmn)
    @inbounds begin
        for k = 1:minmn
            # find index max
            kp = k
            if k < m
                amax = abs(A[k, k])
                for i = k+1:m
                    absi = abs(A[i,k])
                    if absi > amax
                        kp = i
                        amax = absi
                    end
                end
            end
            ipiv[k] = kp
            if !iszero(A[kp,k])
                if k != kp
                    # Interchange
                    for i = 1:n
                        tmp = A[k,i]
                        A[k,i] = A[kp,i]
                        A[kp,i] = tmp
                    end
                end
                # Scale first column
                Akkinv = inv(A[k,k])
                for i = k+1:m
                    A[i,k] *= Akkinv
                end
            elseif info == 0
                info = k
            end
            # Update the rest
            for j = k+1:n
                for i = k+1:m
                    A[i,j] -= A[i,k]*A[k,j]
                end
            end
        end
    end
    #check && checknonsingular(info, pivot)
    return LU(A, ipiv, convert(BlasInt, info))
end

function unsafe_LU!(F::LU)
    return unsafe_LU!(F.factors,F.ipiv)
end

export CholeskyNewton,LUPNewton,static_linsolve
export RestrictedLineSearch

#=
function NLSolvers.find_steplength(mstyle::NLSolvers.MethodStyle, ls::RestrictedLineSearch{F,LS}, φ::T, λ) where {F,LS,T}
    α = ls.f(φ,λ)
    return NLSolvers.find_steplength(mstyle,ls.ls,φ,α)
end =#

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
