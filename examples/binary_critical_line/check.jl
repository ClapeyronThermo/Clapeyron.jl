using Test
include("run.jl")

const constraints = CriticalLineContinuation.NormalizedCriticality.binary_constraints

@testset "Local binary critical-line example" begin
    model = example_model()
    raw = example_trace(model, :raw)
    normalized = example_trace(model)
    @test normalized.bounds == raw.bounds
    for run in (raw, normalized)
        @test length(run.result.branches) == 2
        for (i, (branch, points)) in enumerate(zip(run.result.branches, run.result.points))
            @test branch.status === :coordinate_boundary
            @test length(points) > 1
            endpoint = last(points)
            @test endpoint.T ≈ run.bounds[i] rtol = 2e-10
            @test norm(endpoint.residual, Inf) <= 1e-10
            @test all(isfinite, endpoint.raw_residual)
            @test all(isfinite, endpoint.normalized_residual)
        end
    end
    for (a, b) in zip(last.(raw.result.points), last.(normalized.result.points))
        @test [a.T, a.rho, a.x] ≈ [b.T, b.rho, b.x] rtol = 1e-9
    end

    point = [log(450.0), log(1000.0), 0.5]
    F(u) = constraints(model, exp(-u[2]), exp(u[1]), [u[3], 1 - u[3]]).normalized
    J = Clapeyron.ForwardDiff.jacobian(F, point)
    for h in (1e-4, 1e-5, 1e-6)
        Jfd = hcat(((F(point + h * I(3)[:, i]) - F(point - h * I(3)[:, i])) / (2h)
                   for i in eachindex(point))...)
        @test norm(J - Jfd, Inf) <= 3e-6 * max(1.0, norm(J, Inf))
    end
end
