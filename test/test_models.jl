using OpenSAFT, Test

@testset "ogSAFT" begin
    test_system = system(["methanol"], "ogSAFT")
    Temp        = 298.15
    Vol         = 1e-3
    z           = create_z(test_system, [1.])
    @test typeof(test_system) == OpenSAFT.PCSAFT
    @testset "params" begin
        @test test_system.params.segment[Set(["methanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.sigma[Set(["methanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.epsilon[Set(["methanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.epsilon_assoc[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 1 rtol = 1e10
        @test test_system.params.bond_vol[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 1 rtol = 1e10
        @test test_system.params.n_sites[Set(["methanol"])]["e"] ≈ 1 rtol = 1e10
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_HS(test_system, z, Vol, Temp) ≈ 1 rtol = 1e10
        @test OpenSAFT.a_disp(test_system, z, Vol, Temp) ≈ 1 rtol = 1e10
        @test OpenSAFT.a_chain(test_system, z, Vol, Temp) ≈ 1 rtol = 1e10
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ 1 rtol = 1e10
        @test OpenSAFT.eos(test_system, z, Vol, Temp) ≈ 1 rtol = 1e10
    end
end

@testset "PCSAFT" begin
    test_system = system(["methanol", "ethanol"], "PCSAFT")
    Temp        = 298.15
    Vol         = 1e-3
    z           = create_z(test_system, [0.5, 0.5])
    @test typeof(test_system) == OpenSAFT.PCSAFT
    @testset "params" begin
        @test test_system.params.segment[Set(["ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.sigma[Set(["ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.sigma[Set(["methanol", "ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.epsilon[Set(["ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.epsilon[Set(["methanol", "ethanol"])] ≈ 1 rtol = 1e10
        @test test_system.params.epsilon_assoc[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 1 rtol = 1e10
        @test test_system.params.n_sites[Set(["methanol"])]["e"] ≈ 1 rtol = 1e10
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_HC(test_system, z, Vol, Temp) ≈ 1 rtol = 1e10
        @test OpenSAFT.a_disp(test_system, z, Vol, Temp) ≈ 1 rtol = 1e10
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ 1 rtol = 1e10
        @test OpenSAFT.eos(test_system, z, Vol, Temp) ≈ 1 rtol = 1e10
    end
end
