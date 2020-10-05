using OpenSAFT, Test

@testset "PCSAFT" begin
    test_system = system(["methanol", "ethanol", "cyclohexane"], "PCSAFT")
    @test typeof(test_system) == OpenSAFT.PCSAFT
    @testset "params" begin
        @test test_system.params.segment[Set(["ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.sigma[Set(["ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.sigma[Set(["methanol", "ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.epsilon[Set(["ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.epsilon[Set(["methanol", "ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.epsilon[Set(["methanol", "cyclohexane"])] ≈ 1 rtol = 1e10
        @test test_system.params.epsilon_assoc[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 1 rtol = 1e10
        @test test_system.params.n_sites[Set(["methanol"])]["e"] ≈ 1 rtol = 1e10
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_ideal(test_system, create_z(test_system, [0.4, 0.4, 0.2]), 1, 298) ≈ 1 rtol = 1e10
        @test OpenSAFT.a_res(test_system, create_z(test_system, [0.4, 0.4, 0.2]), 1, 298) ≈ 1 rtol = 1e10
        @test OpenSAFT.eos(test_system, create_z(test_system, [0.4, 0.4, 0.2]), 1, 298) ≈ 1 rtol = 1e10
    end
end

