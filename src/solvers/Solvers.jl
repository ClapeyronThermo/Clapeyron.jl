module Solvers

using LinearAlgebra
using NLSolvers,Roots
using PositiveFactorizations
using DiffResults, ForwardDiff
using StaticArrays

    function solve_cubic_eq(poly::AbstractVector{T}) where {T<:Real}
        # copied from PolynomialRoots.jl, adapted to be AD friendly
        # Cubic equation solver for complex polynomial (degree=3)
        # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
        _1 = one(T)
        third = _1/3
        a1  =  complex(one(T) / poly[4])
        E1  = -poly[3]*a1
        E2  =  poly[2]*a1
        E3  = -poly[1]*a1
        s0  =  E1
        E12 =  E1*E1
        A   =  2*E1*E12 - 9*E1*E2 + 27*E3 # = s1^3 + s2^3
        B   =  E12 - 3*E2                 # = s1 s2
        # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
        Δ = (A*A - 4*B*B*B)^0.5
        if real(conj(A)*Δ)>=0 # scalar product to decide the sign yielding bigger magnitude
            s1 = exp(log(0.5 * (A + Δ)) * third)
        else
            s1 = exp(log(0.5 * (A - Δ)) * third)
        end
        if s1 == 0
            s2 = s1
        else
            s2 = B / s1
        end
        zeta1 = complex(-0.5, sqrt(T(3.0))*0.5)
        zeta2 = conj(zeta1)
        return (third*(s0 + s1 + s2), third*(s0 + s1*zeta2 + s2*zeta1), third*(s0 + s1*zeta1 + s2*zeta2))
    end

    function roots3(x) 
        return SVector(solve_cubic_eq(x))
    end

    function roots3(a,b,c,d) 
        x = SVector(a,b,c,d)
        return roots3(x)
    end


    """
    function x_sol(res::NLSolvers.ConvergenceInfo)
        
        Returns the scalar or vector x that solves the system of equations or is the minimizer of an optimization procedure.
    """

    function x_sol(res::NLSolvers.ConvergenceInfo{<:NLSolvers.LineSearch, <:Any, <:NLSolvers.NEqOptions})
        res.info.solution
    end

    function x_sol(res::NLSolvers.ConvergenceInfo{<:NLSolvers.TrustRegion, <:Any, <:NLSolvers.NEqOptions})
        return res.info.zero
    end

    function x_sol(res::NLSolvers.ConvergenceInfo{<:Any, <:Any, <:NLSolvers.OptimizationOptions})
        return res.info.minimizer
    end

   
    include("ADNewton.jl")
    include("nested.jl")
    include("nlsolve.jl")
    include("fixpoint/fixpoint.jl")
    include("optimize.jl")
    #include("tunneling.jl") work in progress
end # module
