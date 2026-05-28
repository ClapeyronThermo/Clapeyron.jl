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
