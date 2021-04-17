module Solvers

using LinearAlgebra
using NLopt, NLSolvers,Roots
#using NLSolve
using  DiffResults, ForwardDiff
using StaticArrays
using PolynomialRoots
include("tunneling.jl")
include("ADNewton.jl")
include("nested.jl")
include("nlsolve.jl")
include("fixpoint/fixpoint.jl")

polyroots(x) = PolynomialRoots.roots(x,polish=true)

end # module
