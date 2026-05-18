module ClapeyronMetaheuristicsExt

using Clapeyron
using Metaheuristics

function logger(st)
    if st.iteration % 10 == 0
        println("Iteration: ",st.iteration)
        println("Best solution: ",st.best_sol.x)
        println("Best value: ",st.best_sol.f)
    end
end

function Metaheuristics.optimize(f::Function,estimator::Estimation{T},method::Metaheuristics.AbstractAlgorithm = Metaheuristics.ECA();verbose=false,logger::Function = (status) -> logger(status)) where T<:EoSModel
    nparams = length(estimator.toestimate.upper)
    guesses = [estimator.toestimate.guess[i][1] for i ∈ 1:nparams]
    Metaheuristics.set_user_solutions!(method, guesses, f)
    upper = [estimator.toestimate.upper[i][1] for i ∈ 1:nparams]
    lower = [estimator.toestimate.lower[i][1] for i ∈ 1:nparams]
    bounds  = [lower upper]'
    if verbose
        results = Metaheuristics.optimize(f,bounds,method;logger=logger)
    else
        results = Metaheuristics.optimize(f,bounds,method;logger=(status) -> nothing)
    end
    
    params = minimizer(results)
    model = return_model(estimator,estimator.model,params)
    return params, model
end

end