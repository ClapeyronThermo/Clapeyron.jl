using LinearAlgebra
using Simplex

#A = Float64[1 -1 1 1 0 0; 1 1 0 0 -1 0; 0 0 1 0 0 1]
#b = Float64[4, 0, 6]
#c = Float64[1, 2, 3, 0, 0, 0]

A =  Float64[3 2 1 2 1 0 0; 1 1 1 1 0 1 0; 4 3 3 4 0 0 1]
b =  Float64[225, 117, 420]
c = -Float64[19, 13, 12, 17, 0, 0, 0]

(primal_sol, dual_sol, obj) = fullrsm(A, b, c)

println("primal sol         = $(primal_sol)")
println("dual sol           = $(dual_sol)")
println("objective function = $(-obj)")


#include("/home/cbranch/julia/dev/simplex.jl/simplex.jl")

#simplex()
