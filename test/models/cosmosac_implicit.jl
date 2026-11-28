@testset "COSMOSAC implicit differentiation" begin
    FD = Clapeyron.ForwardDiff
    # One segment has the exact solution Γ = PW^(-1/2). Nested Duals must
    # include the indirect lower-order terms, not just a final Newton step.
    root(x) = Clapeyron.cosmosac_solve(reshape([x], 1, 1), 1000, 1e-14)[1]
    x = 2.3
    @test root(x) ≈ x^(-1/2) rtol=1e-13
    @test FD.derivative(root, x) ≈ -x^(-3/2)/2 rtol=1e-12
    @test FD.derivative(t -> FD.derivative(root, t), x) ≈ 3x^(-5/2)/4 rtol=1e-11
    @test FD.derivative(t -> FD.derivative(u -> FD.derivative(root, u), t), x) ≈ -15x^(-7/2)/8 rtol=1e-10

    for Model in (COSMOSAC02, COSMOSAC10, COSMOSACdsp)
        model = Model(["water", "ethanol"], puremodel=BasicIdeal())
        # Scale temperature so finite-difference directions have similar size.
        f(q) = Clapeyron.excess_gibbs_free_energy(model, 1e5, 300q[1], q[2:end])
        q = [333.15/300, 0.4, 0.6]
        gradient = FD.gradient(f, q)
        hessian = FD.hessian(f, q)
        h = 1e-3
        for i in eachindex(q)
            ei = zeros(length(q))
            ei[i] = h
            gfd = (f(q-2ei)-8f(q-ei)+8f(q+ei)-f(q+2ei))/(12h)
            hfd = (-f(q+2ei)+16f(q+ei)-30f(q)+16f(q-ei)-f(q-2ei))/(12h^2)
            @test gradient[i] ≈ gfd rtol=2e-6 atol=2e-6
            @test hessian[i, i] ≈ hfd rtol=2e-5 atol=2e-5
            for j in 1:(i - 1)
                ej = zeros(length(q))
                ej[j] = h
                mixed = (f(q+ei+ej)-f(q+ei-ej)-f(q-ei+ej)+f(q-ei-ej))/(4h^2)
                @test hessian[i, j] ≈ mixed rtol=1e-4 atol=2e-4
            end
        end
        @test hessian ≈ transpose(hessian) rtol=1e-10 atol=1e-9
        # Fixed-support differentiation at pure compositions remains valid in T.
        @test isfinite(FD.derivative(T -> Clapeyron.excess_gibbs_free_energy(model, 1e5, T, [0.0, 1.0]), 333.15))
    end

    # Exercise support selection itself: water contributes bins absent in benzene.
    model = COSMOSAC02(["water", "benzene"], puremodel=BasicIdeal())
    A = model.params.A.values
    Pi = model.params.Pi.values
    profile(x) = (x*Pi[1]+(1-x)*Pi[2])/(x*A[1]+(1-x)*A[2])
    P = profile(0.0)
    Pd = profile(FD.Dual(0.0, 1.0))
    primal_support = findall(!iszero, P)
    dual_support = findall(!iszero, Pd)
    @test length(dual_support) > length(primal_support)
    PW = Clapeyron.PΔW(model, 1e5, 333.15, [0.0, 1.0], Pd, dual_support)
    Γ = Clapeyron.cosmosac_solve(Clapeyron.nested_pvalue(PW), 1000, 1e-12)
    lnΓ = Clapeyron._lnΓ(model, 1e5, 333.15, [FD.Dual(0.0, 1.0), FD.Dual(1.0, -1.0)], Pd)
    @test exp.(Clapeyron.nested_pvalue.(lnΓ[dual_support])) ≈ Γ rtol=1e-12
    @test all(>(0), Γ)

    # A zero primal column with nonzero partials must survive in the reduced
    # matrix. Its segment coefficient is needed when the component appears.
    W(t) = [2.0 t; 1.0 3t]
    roots(t) = Clapeyron.cosmosac_solve(W(t), 1000, 1e-14)
    Γ = roots(0.0)
    @test Γ ≈ [1/sqrt(2), sqrt(2)] rtol=1e-13
    J = [2sqrt(2) 0.0; sqrt(2) 1/sqrt(2)]
    expected = -(J \ (Γ .* ([0.0 1.0; 0.0 3.0]*Γ)))
    @test FD.jacobian(t -> roots(t[1]), [0.0])[:, 1] ≈ expected rtol=1e-11
end
