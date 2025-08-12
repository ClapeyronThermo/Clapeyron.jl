module TRS

using LinearAlgebra
using LinearMaps
using IterativeSolvers
using Polynomials

include("eigenproblems.jl")
include("trust_region_boundary.jl")
include("trust_region_small.jl")
export trs_small, trs_boundary_small

end
