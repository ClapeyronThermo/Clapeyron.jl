using Clapeyron, LinearAlgebra
include("CriticalLineContinuation.jl")
using .CriticalLineContinuation: ContinuationOptions, trace_binary_critical_line

function example_model()
    return PCSAFT(["methane", "decane"]; userlocations = (;
        segment = [1.0, 4.6627], sigma = [3.7039, 3.8384],
        epsilon = [150.03, 243.87], Mw = [16.04, 142.29],
        k = zeros(2, 2), l = zeros(2, 2),
    ))
end

function example_trace(model, residual_mode = :normalized)
    # Retained model-equation root, not an experimental observation or fit.
    x = 0.5
    T = 579.2905095417091
    V = inv(3029.3162302379783)
    p = pressure(model, V, T, [x, 1 - x])
    bounds = (0.999 * T, 1.001 * T)
    result = trace_binary_critical_line(
        model; x0 = x, seed = (T, p, V), Tbounds = bounds,
        stop_coordinate = :temperature, options = ContinuationOptions(max_steps = 5),
        residual_mode,
    )
    return (; result, bounds)
end

function main(io = stdout)
    model = example_model()
    result = example_trace(model).result
    for (branch, points) in zip(result.branches, result.points)
        endpoint = last(points)
        println(io, (; status = branch.status, endpoint.T, endpoint.x, endpoint.p,
            normalized_residual = norm(endpoint.normalized_residual, Inf)))
    end
    return nothing
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
