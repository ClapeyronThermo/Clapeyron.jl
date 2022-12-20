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
    _nlsolve(f!,x0,method,options,chunk)
end

function _nlsolve(f!,x0,method,options::NEqOptions,chunk)
    vector_objective = autoVectorObjective(f!,x0,chunk)
    nl_problem = NEqProblem(vector_objective)
    return nlsolve(nl_problem, x0,method, options)
end

function _nlsolve(f!::F,x0,method::NLSolvers.Anderson,options::NEqOptions,chunk) where F
    function g!(out,x)
        f!(out,x)
        out .+= x
        return out
    end
    __fixpoint(g!,x0,method,options.f_abstol,0.0,options.maxiter,true)
end

function nlsolve(nl_problem::NEqProblem,x0,method =TrustRegion(Newton(), NWI()),options=NEqOptions())
    return _nlsolve(nl_problem, x0, method, options)
end

function _nlsolve(nl_problem::NEqProblem,x0,method,options::NEqOptions)
    return NLSolvers.solve(nl_problem, x0, method, options)
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
x_converged(ci::NLSolvers.ConvergenceInfo) = false
f_converged(ci::NLSolvers.ConvergenceInfo) = false
g_converged(ci::NLSolvers.ConvergenceInfo) = false

function x_converged(ci::NLSolvers.ConvergenceInfo{<:NLSolvers.NelderMead,<:Any,<:NLSolvers.OptimizationOptions})
    return ci.info.ρs <= ci.options.x_abstol
end

function x_converged(ci::NLSolvers.ConvergenceInfo{<:NLSolvers.SimulatedAnnealing,<:Any,<:NLSolvers.OptimizationOptions})
    #don't know what to do do here, simulated annealing does not "converge" in the typical sense
    return true
end

function x_converged(ci::NLSolvers.ConvergenceInfo{<:Any,<:Any,<:NLSolvers.OptimizationOptions})
    info = ci.info
    opt = ci.options
    x_abstol = haskey(info, :ρs) && (info.ρ <= opt.x_abstol)
    x_reltol = haskey(info, :ρs) && (info.ρs/info.ρx <= opt.x_reltol)
    return x_abstol || x_reltol
end

function f_converged(ci::NLSolvers.ConvergenceInfo{<:Any,<:Any,<:NLSolvers.OptimizationOptions})
    info = ci.info
    opt = ci.options
    f_limit = !isfinite(opt.f_limit) || (info.minimum <= opt.f_limit)
    f_abstol = haskey(info, :fx) && (abs(info.fx - info.minimum) <= opt.f_abstol)
    f_reltol = haskey(info, :fx) && (abs((info.fx - info.minimum)/info.fx) <= opt.f_reltol)
    return f_limit || f_abstol || f_reltol
end

function f_converged(ci::NLSolvers.ConvergenceInfo{<:Any,<:Any,<:NLSolvers.NEqOptions})
    info = ci.info
    opt = ci.options
    ρF = norm(ci.info.best_residual, Inf)
    #ρFz = norm(ci.info.solution, 2)
    f_abstol = ρF <= opt.f_abstol
    return f_abstol
end

function g_converged(ci::NLSolvers.ConvergenceInfo{<:Any,<:Any,<:NLSolvers.OptimizationOptions})
    info = ci.info
    opt = ci.options
    if haskey(info, :∇fz)
        ρ∇f = opt.g_norm(info.∇fz)
        g_abstol = ρP<=opt.g_abstol
        g_reltol = ρ∇f/info.∇f0<=opt.g_reltol
        if haskey(info, :prob) && hasbounds(info.prob)
            ρP = opt.g_norm(
                info.solution .- clamp.(info.solution .- info.∇fz, info.prob.bounds...),
            )
            gp_abstol = ρP <= opt.g_abstol
        else
            gp_abstol = false
        end
    else
        g_abstol =  false
        g_reltol = false
    end
    return g_abstol || g_reltol || gp_abstol
end


function converged(ci::NLSolvers.ConvergenceInfo)
    opt = ci.options
    info = ci.info
    conv_flags = x_converged(ci) || f_converged(ci) || g_converged(ci)
    if haskey(info, :Δ)
        Δmin = ci.solver.Δupdate.Δmin isa Nothing ? 0 : ci.solver.Δupdate.Δmin
        Δ = info.Δ <= Δmin
        Δ_finite = isfinite(info.Δ)
    else
        Δ = false
        Δ_finite = true
    end
    x = x_converged(ci)
    f = f_converged(ci)
    g = g_converged(ci)
    #finite flags
    x_finite = all(isfinite,x_sol(ci))
    conv_flags = f || x || g || Δ
    finite_flags = x_finite && Δ_finite #  && g_finite && f_finite
    return conv_flags && finite_flags
end 
