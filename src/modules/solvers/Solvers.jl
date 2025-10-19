module Solvers

using LinearAlgebra
using LinearAlgebra: BlasInt
using NLSolvers,Roots
using PositiveFactorizations
using DiffResults, ForwardDiff
using StaticArrays
using Roots


export CholeskyNewton,static_linsolve,Newton2

__is_implace(x::Number) = false
__is_implace(x::Array) = true
__is_implace(x::MVector) = true
__is_implace(x::SVector) = false
__is_implace(x::StaticArrays.SizedVector) = true

#__is_implace(x::AbstractVector) = ArrayInterface.can_setindex(x)

"""
    solution(res::NLSolvers.ConvergenceInfo)
    
Returns the scalar or vector x that solves the system of equations or is the minimizer of an optimization procedure.
"""
solution(res) = NLSolvers.solution(res)
function solution(res::NLSolvers.ConvergenceInfo{NLSolvers.BrentMin{Float64}})
    return res.info.x
end

const x_sol = solution

x_minimum(res::NLSolvers.ConvergenceInfo) = res.info.minimum
#for BrentMin (should be fixed at NLSolvers 0.3)
x_minimum(res::Tuple{<:Number,<:Number}) = last(res)


include("linsolve.jl")
include("nlsolvers.jl")
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
include("npipm.jl")

end # module
