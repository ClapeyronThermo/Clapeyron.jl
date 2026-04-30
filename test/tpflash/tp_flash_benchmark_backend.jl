"""
Optimizer backend for TP flash benchmark.

This indirection exists because (per project workflow) switching metaheuristics is
done by hacking source code. Edit this file to change the optimizer used by
`run_case_with_optimizer` without touching the benchmark logic.
"""

using Clapeyron
try
    using BlackBoxOptim
catch
    # Optional dependency: only needed when `backend == :bbo`.
end

"""
Backend optimization hook for the TP-flash benchmark.

Returns `(best_u, best_y)` for the provided scalar objective.

Supported backends:
- `:rdex`  -> `Clapeyron.Solvers.RDEx` (ask–tell, supports stagnation early stop)
- `:sass`  -> `Clapeyron.Solvers.SASS` (ask–tell, supports stagnation early stop)
- `:bbo`   -> `BlackBoxOptim.bboptimize` (batch optimizer, uses `MaxFuncEvals`)

This indirection exists because (per project workflow) switching metaheuristics is
done by hacking source code. Edit this file if you want to change the mapping
from `backend` to a concrete optimizer.
"""
function tpflash_optimize(
    objective::Function;
    backend::Symbol,
    population_size::Int,
    fe_max::Int,
    lb::Vector{Float64},
    ub::Vector{Float64},
    seed::Int,
    time_limit::Float64,
    stagnation_evals::Int,
    stagnation_tol::Float64,
    verbose::Bool=false,
)
    if backend == :sass
        algo = Clapeyron.Solvers.SASS(population_size, fe_max, lb, ub; seed, stagnation_evals, stagnation_tol)
        deadline_ns = isfinite(time_limit) ? time_ns() + floor(Int, 1e9time_limit) : typemax(Int)
        while !Clapeyron.Solvers.isdone(algo) && algo.state.nfes < fe_max
            u = Clapeyron.Solvers.ask!(algo)
            y = objective(u)
            Clapeyron.Solvers.tell!(algo, y)
            time_ns() >= deadline_ns && break
        end
        best_u, best_y = Clapeyron.Solvers.best(algo)
        return collect(Float64, best_u), Float64(best_y)
    elseif backend == :bbo
        isdefined(@__MODULE__, :BlackBoxOptim) || error("backend :bbo requires BlackBoxOptim; it is not available in this environment")
        max_time = isfinite(time_limit) ? Float64(time_limit) : 0.0
        ranges = [(lb[i], ub[i]) for i in eachindex(lb)]

        res = BlackBoxOptim.bboptimize(
            objective;
            SearchRange=ranges,
            NumDimensions=length(lb),
            MaxFuncEvals=fe_max,
            MaxTime=max_time,
            PopulationSize=population_size,
            RandomSeed=seed,
            TraceMode=verbose ? :verbose : :silent,
            # Best-effort mapping of stagnation controls (not strictly equivalent to RDEx):
            MaxStepsWithoutProgress=stagnation_evals > 0 ? stagnation_evals : 0,
            MinDeltaFitnessTolerance=stagnation_tol,
        )
        return BlackBoxOptim.best_candidate(res), Float64(BlackBoxOptim.best_fitness(res))
    else
        error("unknown tpflash backend: $backend (expected :rdex, :sass, or :bbo)")
    end
end
