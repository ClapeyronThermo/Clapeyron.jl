#nlsolve functionality
"""
    function nlsolve(f!,x0,method=TrustRegion(Newton(), NWI()), options=NEqOptions(),chunk = ForwardDiff.Chunk{2}())


Given a function `f!(result,x)` that returns a system of equations,
`nlsolve(f!,x0)` returns a `NLSolvers.ConvergenceInfo` struct that contains the results of the non-linear solving procedure.

Uses `NLSolvers.jl` as backend, the jacobian is calculated with `ForwardDiff.jl`, with the specified `chunk` size

To obtain the underlying solution vector, use [`x_sol`](@ref)

To see available solvers and options, check `NLSolvers.jl`
"""
function nlsolve(f!,x0,method=TrustRegion(Newton(), NWI()),options=NEqOptions(),chunk = ForwardDiff.Chunk{2}())
    if method isa Roots.AbstractUnivariateZeroMethod
        return roots_nlsolve(f!,x0,method,options)
    end
    vector_objective = autoVectorObjective(f!,x0,chunk)
    nl_problem = NEqProblem(vector_objective; inplace = _inplace(x0))
    return nlsolve(nl_problem, x0,method, options)
end

function nlsolve(nl_problem::NEqProblem,x0,method =TrustRegion(Newton(), NWI()),options=NEqOptions())
    return NLSolvers.solve(nl_problem, x0,method, options)
end

#Roots.jl support
function roots_nlsolve(f::F,x0,method::Roots.AbstractBracketingMethod,options) where F
    prob = Roots.ZeroProblem(f,x0)
    brk = f(x0[1])*f(x0[2])
    if brk > 0
        sol = zero(brk)/zero(brk)
    else
        sol = Roots.solve(prob,method)
    end    
end

function roots_nlsolve(f::F,x0::Number,method::Roots.AbstractNonBracketingMethod ,options) where F
    prob = Roots.ZeroProblem(f,x0)
    sol = Roots.solve(prob,method)
end

function roots_nlsolve(f::F,x0::Number,method::Roots.AbstractNewtonLikeMethod ,options) where F
    function fdf(vz)
        fx,dfx = Solvers.f∂f(f,vz)
        return fx,fx/dfx
    end
    prob = Roots.ZeroProblem(fdf,x0)
    sol = Roots.solve(prob,method)
end

function roots_nlsolve(f::F,x0::Number,method::Roots.AbstractHalleyLikeMethod,options) where F
    function d2f(vz)
        fx,dfx,d2fx = Solvers.f∂f∂2f(f,vz)
        return fx,fx/dfx,dfx/d2fx
    end
    prob = Roots.ZeroProblem(d2f,x0)
    sol = Roots.solve(prob,method)
end


function autoVectorObjective(f!,x0,chunk)
    Fcache = x0 .* false
    jconfig = ForwardDiff.JacobianConfig(f!,x0,x0,chunk)
    function j!(J,x)
        ForwardDiff.jacobian!(J,f!,Fcache,x,jconfig)
        J
    end
    function fj!(F,J,x)
        ForwardDiff.jacobian!(J,f!,F,x,jconfig)
        F,J
    end
    return NLSolvers.VectorObjective(f!,j!,fj!,nothing)
end

_inplace(x0) = true
_inplace(x0::SVector) = false

function autoVectorObjective(f!,x0::StaticArrays.SVector{2,T},chunk) where T
    f(x) = f!(nothing,x) #we assume that the F argument is unused in static arrays
    j(J,x) = ForwardDiff.jacobian(f,x)
    fj(F,J,x) = FJ_ad(f,x)
    return NLSolvers.VectorObjective(f!,j,fj,nothing)
end

function autoVectorObjective(f!,x0::StaticArrays.SVector{3,T},chunk) where T
    f(x) = f!(nothing,x) #we assume that the F argument is unused in static arrays
    j(J,x) = ForwardDiff.jacobian(f,x)
    fj(F,J,x) = FJ_ad(f,x)
    return NLSolvers.VectorObjective(f!,j,fj,nothing)
end

function autoVectorObjective(f!,x0::StaticArrays.SVector,chunk)
    f(x) = f!(nothing,x) #we assume that the F argument is unused in static arrays
    j(J,x) = ForwardDiff.jacobian(f,x)
    fj(F,J,x) = FJ_ad(f,x)
    return NLSolvers.VectorObjective(f,j,fj,nothing)
end

#= only_fj!: NLsolve.jl legacy form:

function only_fj!(F, J, x)
    # shared calculations begin
    # ...
    # shared calculation end
    if !(F == nothing)
        # mutating calculations specific to f! goes here
    end
    if !(J == nothing)
        # mutating calculations specific to j! goes
    end
end
=#
function only_fj!(fj!::T) where T
    function _f!(F,x)
        fj!(F,nothing,x)
        F
    end

    function _fj!(F,J,x)
        fj!(F,J,x)
        F,J
    end

    function _j!(J,x)
        fj!(nothing,J,x)
        J
    end
    return NLSolvers.VectorObjective(_f!,_j!,_fj!,nothing) |> NEqProblem
    # return NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> NEqProblem
end

#=
#trying to make nlsolve(f,x0,LineSearch(Newton(),HZAW())) work

function NLSolvers.upto_gradient(meritobj::NLSolvers.MeritObjective, ∇f, x)
    neq = meritobj.prob
    G = neq.R.F(∇f, x)
    F = (norm(G)^2) / 2
    return F,G
end
=#

#=

=#

struct Newton2Var end

function nlsolve2(f::FF,x::SVector{NN,TT},method::Newton2Var,options=NEqOptions()) where {FF,NN,TT}
    function FJ(_z)
        return FJ_ad(f,_z)
    end
    Fx, Jx = FJ(x)
    z = x
    T = eltype(Fx)
    stoptol = T(options.f_abstol)
    ρF0, ρ2F0 = norm(Fx, Inf), norm(Fx, 2)
    nan = T(NaN)
    ρs = nan
    #@show ρF0
    if ρF0 < stoptol
        return x
    end
    iter = 1
    converged = false
    while iter ≤ options.maxiter
        d = Jx \ -Fx
        #@show Jx, Fx
        x = x + d
        Fx, Jx = FJ(x)
        ρF = norm(Fx, Inf)
        ρs = norm(d, Inf)
        #@show ρF, ρs
        if ρs <= stoptol || ρF <= stoptol
            converged = true
            break
        end

        if !all(isfinite,x)
            converged = false
            break
        end
        iter += 1
    end
    if !converged
        x  = nan .* x
    end
    return x
end
