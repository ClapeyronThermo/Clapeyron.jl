
@testset "SAFT models" begin
    T = 298.15
    V = 1e-4

    @testset "ogSAFT" begin
        system = ogSAFT(["water","ethylene glycol"])
        z = [0.5, 0.5]
        @test Clapeyron.a_seg(system, V, T, z) ≈ -2.0332062924093366 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ -0.006317441684202759 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -4.034042081699316 rtol = 1e-6
    end

    @testset "CKSAFT" begin
        system = CKSAFT(["carbon dioxide", "2-propanol"])
        z = [0.5, 0.5]
        @test Clapeyron.a_seg(system, V, T, z) ≈ -1.24586302917188 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ -0.774758615408493 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.2937079004096872 rtol = 1e-6
    end

    @testset "sCKSAFT" begin
        system = sCKSAFT(["benzene","acetic acid"])
        z = [0.5, 0.5]
        @test Clapeyron.a_seg(system, V, T, z) ≈ -3.1809330810925256 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -3.3017434376105514 rtol = 1e-6
    end

    @testset "BACKSAFT" begin
        system = BACKSAFT(["carbon dioxide"])
        z = [1.]
        @test Clapeyron.a_hcb(system, V, T, z) ≈ 1.0118842111801198 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ -0.14177009317268635 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z) ≈ -2.4492518566426296 rtol = 1e-6
    end

    @testset "CPA" begin
        system = CPA(["ethanol","benzene"])
        z = [0.5, 0.5]
        @test Clapeyron.a_SRK(system, V, T, z) ≈ 4.510022402195623 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.1575210505284332 rtol = 1e-6
    end

    @testset "SAFTVRSW" begin
        system = SAFTVRSW(["water", "ethane"])
        z = [0.5, 0.5]
        @test Clapeyron.a_mono(system, V, T, z) ≈ -1.4367205951569462 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ 0.024000058201261557 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -0.5238154638538838 rtol = 1e-6
    end

    @testset "softSAFT" begin
        system = softSAFT(["hexane","1-propanol"])
        z = [0.5,0.5]
        @test Clapeyron.a_LJ(system, V, T, z) ≈ -3.960728242264164 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ 0.3736728407455211 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -2.0461376618069034 rtol = 1e-6
    end

    @testset "PCSAFT" begin
        system = PCSAFT(["butane", "ethanol"])
        z = [0.5, 0.5]
        @test Clapeyron.a_hc(system, V, T, z) ≈ 3.1148229872928654 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z) ≈ -6.090736508783152 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.6216064387201956 rtol = 1e-6
    end

    @testset "sPCSAFT" begin
        system = sPCSAFT(["propane", "methanol"])
        z = [0.5, 0.5]
        @test Clapeyron.a_hc(system, V, T, z) ≈ 2.024250583187793 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z) ≈ -4.138653131750594 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.1459701721909195 rtol = 1e-6
        #difference in this error is almost exactly 0.5, suspicious
    end

    @testset "SAFTVRMie" begin
        system = SAFTVRMie(["methanol", "water"])
        z = [0.5, 0.5]
        @test Clapeyron.a_mono(system, V, T, z) ≈ -0.9729134860869052 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ -0.02834738013535014 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -4.180807072390184 rtol = 1e-6
    end

    @testset "SAFTVRQMie" begin
        system = SAFTVRQMie(["helium"])
        z = [1.]
        @test Clapeyron.a_mono(system, V, T, z) ≈ 0.12286776703976324 rtol = 1e-6
    end

    @testset "SAFTgammaMie" begin
        system = SAFTgammaMie(["ethanol"])
        V_γMie = exp10(-3.5)
        z = [1.]
        @test Clapeyron.a_mono(system, V_γMie, T, z) ≈ -1.151043781769667 rtol = 1e-6
        @test Clapeyron.a_chain(system, V_γMie, T, z) ≈ -0.1255227354789658 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V_γMie, T, z) ≈ -1.9386416653191778 rtol = 1e-6
    end
end