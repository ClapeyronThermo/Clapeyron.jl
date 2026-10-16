module ClapeyronMetaheuristicsExt

using Clapeyron
using Metaheuristics

function logger(st,verbose = false)
    if verbose
        if st.iteration % 10 == 0
            println("Iteration: ",st.iteration)
            println("Best solution: ",st.best_sol.x)
            println("Best value: ",st.best_sol.f)
        end
    end
end

function Metaheuristics.optimize(f::Function,estimator::EstimationProblem,method::Metaheuristics.AbstractAlgorithm = Metaheuristics.ECA();verbose=false,logger = Base.Fix2(logger,verbose))
    guesses = Clapeyron.EstimationUtils.initial_guess(estimator)
    Metaheuristics.set_user_solutions!(method, guesses, f)
    ub = Clapeyron.EstimationUtils.upper_bounds(estimator)
    lb = Clapeyron.EstimationUtils.lower_bounds(estimator)
    bounds = Metaheuristics.boxconstraints(lb = lb, ub = ub)
    results = Metaheuristics.optimize(f,bounds,method;logger=logger)
    params = minimizer(results)
    model = return_model(estimator,estimator.model,params)
    return params, model
end

function Metaheuristics.optimize(estimator::EstimationProblem,method::Metaheuristics.AbstractAlgorithm = Metaheuristics.ECA();verbose=false,logger = Returns(nothing))
    f = Clapeyron.EstimationUtils.objective_function(estimator)
    return Metaheuristics.optimize(f,estimator,method;verbose,logger)
end

end #module