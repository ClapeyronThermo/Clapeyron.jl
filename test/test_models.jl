
@testset "models" begin
    T = 298.15
    V = 1e-4

    @testset "ogSAFT" begin
        system = ogSAFT(["methanol"])
        z = [1.]
        @test OpenSAFT.a_seg(system, V, T, z) ≈ -0.540464960274975 rtol = 1e-6
        @test OpenSAFT.a_chain(system, V, T, z) ≈ -0.22445414690633522 rtol = 1e-6
        @test OpenSAFT.a_assoc(system, V, T, z) ≈ -4.695916382523908 rtol = 1e-6
    end

    @testset "CKSAFT" begin
        system = CKSAFT(["carbon dioxide", "2-propanol"])
        z = [0.5, 0.5]
        @test OpenSAFT.a_seg(system, V, T, z) ≈ -1.2395529662948277 rtol = 1e-6
        @test OpenSAFT.a_chain(system, V, T, z) ≈ -0.7747586154084931 rtol = 1e-6
        @test OpenSAFT.a_assoc(system, V, T, z) ≈ -1.7937079004096872 rtol = 1e-6
    end

    @testset "SAFTVRSW" begin
        system = SAFTVRSW(["water", "ethane"])
        z = [0.5, 0.5]
        @test OpenSAFT.a_mono(system, V, T, z) ≈ -1.4659622048306407 rtol = 1e-6
        @test OpenSAFT.a_chain(system, V, T, z) ≈ 0.022703334973543182 rtol = 1e-6
        @test OpenSAFT.a_assoc(system, V, T, z) ≈ -0.5091186885233859 rtol = 1e-6
    end

    @testset "softSAFT" begin
        system = softSAFT(["methanol"])
        z = [1.]
        @test OpenSAFT.a_LJ(system, V, T, z) ≈ -1.299697047509115 rtol = 1e-6
        @test OpenSAFT.a_chain(system, V, T, z) ≈ 0.33553325545339646 rtol = 1e-6
        @test OpenSAFT.a_assoc(system, V, T, z) ≈ -3.7490459421490447 rtol = 1e-6
    end

    @testset "PCSAFT" begin
        system = PCSAFT(["butane", "ethanol"])
        z = [0.5, 0.5]
        @test OpenSAFT.a_hc(system, V, T, z) ≈ 3.114823074155765 rtol = 1e-6
        @test OpenSAFT.a_disp(system, V, T, z) ≈ -6.090736624622517 rtol = 1e-6
        @test OpenSAFT.a_assoc(system, V, T, z) ≈ -2.121606453473655 rtol = 1e-6
    end

    @testset "sPCSAFT" begin
        system = sPCSAFT(["pentane", "methanol"])
        z = [0.5, 0.5]
        @test OpenSAFT.a_hc(system, V, T, z) ≈ 3.568650770403549 rtol = 1e-6
        @test OpenSAFT.a_disp(system, V, T, z) ≈ -6.994181358803752 rtol = 1e-6
        @test OpenSAFT.a_assoc(system, V, T, z) ≈ -1.7525112985184315 rtol = 1e-6
    end

    @testset "SAFTVRMie" begin
        system = SAFTVRMie(["methanol"])
        z = [1.]
        @test OpenSAFT.a_mono(system, V, T, z) ≈ -1.7176380421592838 rtol = 1e-6
        @test OpenSAFT.a_chain(system, V, T, z) ≈ -0.030259270092795967 rtol = 1e-6
        @test OpenSAFT.a_assoc(system, V, T, z) ≈ -3.1565551121889293 rtol = 1e-6
    end

    @testset "SAFTVRQMie" begin
        system = SAFTVRQMie(["helium"])
        z = [1.]
        @test OpenSAFT.a_mono(system, V, T, z) ≈ 0.12253715358076675 rtol = 1e-6
    end

    @testset "SAFTgammaMie" begin
        system = SAFTgammaMie(["ethanol"])
        V_γMie = exp10(-3.5)
        z = [1.]
        @test OpenSAFT.a_mono(system, V_γMie, T, z) ≈ -1.151043781769667 rtol = 1e-6
        @test OpenSAFT.a_chain(system, V_γMie, T, z) ≈ -0.1255227354789658 rtol = 1e-6
        @test OpenSAFT.a_assoc(system, V_γMie, T, z) ≈ -1.9386416653191778 rtol = 1e-6
    end
end
