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
    vector_objective = autoVectorObjective(f!,x0,chunk)
    nl_problem = NEqProblem(vector_objective)
    return nlsolve(nl_problem, x0,method, options)
end

function nlsolve(nl_problem::NEqProblem,x0,method =TrustRegion(Newton(), NWI()),options=NEqOptions())
    return NLSolvers.solve(nl_problem, x0,method, options)
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
    function jv!(x)
        return nothing
    end
    return NLSolvers.VectorObjective(f!,j!,fj!,jv!)
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

    _jv!(x) = nothing
    return NLSolvers.VectorObjective(_f!,_j!,_fj!,_jv!) |> NEqProblem
    # return NLSolvers.VectorObjective(f!,j!,fj!,jv!) |> NEqProblem
end

#= 
#trying to make nlsolve(f,x0,LineSearch(Newton(),HZAW())) work

function NLSolvers.upto_gradient(meritobj::NLSolvers.MeritObjective, ∇f, x)
    neq = meritobj.prob
    G = neq.R.F(∇f, x)
    F =  (norm(G)^2) / 2
    return F,G
end
=#

