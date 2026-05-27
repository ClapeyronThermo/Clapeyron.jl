struct EstimationSolution{L,X,M}
    loss::L
    x::X
    model::M
end

module EstimationUtils
    function parameter_vector end
    function lower_bounds end
    function upper_bounds end
    function initial_guess end
    function objective_function end

    export parameter_vector
    export lower_bounds
    export upper_bounds
    export initial_guess
    export objective_function
end

EstimationUtils.lower_bounds(model::ToEstimate) = model.lower
EstimationUtils.upper_bounds(model::ToEstimate) = model.upper
EstimationUtils.initial_guess(model::ToEstimate) = model.guess

EstimationUtils.lower_bounds(model::EstimationModel) = EstimationUtils.lower_bounds(model.toestimate)
EstimationUtils.upper_bounds(model::EstimationModel) = EstimationUtils.upper_bounds(model.toestimate)
EstimationUtils.initial_guess(model::EstimationModel) = EstimationUtils.initial_guess(model.toestimate)

EstimationUtils.lower_bounds(model::Estimation) = EstimationUtils.lower_bounds(model.toestimate)
EstimationUtils.upper_bounds(model::Estimation) = EstimationUtils.upper_bounds(model.toestimate)
EstimationUtils.initial_guess(model::Estimation) = EstimationUtils.initial_guess(model.toestimate)

EstimationUtils.parameter_vector(model::Estimation) = EstimationUtils.parameter_vector(model.toestimate)
EstimationUtils.parameter_vector(est_model::EstimationModel) = get_eos_parameters(est_model)
