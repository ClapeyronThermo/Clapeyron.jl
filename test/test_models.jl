using OpenSAFT, Test

@testset "ogSAFT" begin
    test_system = system(["methanol"], "ogSAFT")
    Temp        = 298.15
    Vol         = 1e-4
    z           = create_z(test_system, [1.])
    @test typeof(test_system) == OpenSAFT.ogSAFT
    @testset "params" begin
        @test test_system.params.segment[Set(["methanol"])] ≈ 1.6 rtol = 1e-8
        @test test_system.params.sigma[Set(["methanol"])] ≈ 3.203e-10 rtol = 1e-8
        @test test_system.params.epsilon[Set(["methanol"])] ≈ 163.15 rtol = 1e-8
        @test test_system.params.epsilon_assoc[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 2964 rtol = 1e-8
        @test test_system.params.bond_vol[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 0.053 rtol = 1e-8
        @test test_system.params.n_sites[Set(["methanol"])]["e"] ≈ 1 rtol = 1e-8
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_seg(test_system, z, Vol, Temp) ≈ -0.540464960274975 rtol = 1e-8
        @test OpenSAFT.a_chain(test_system, z, Vol, Temp) ≈ -0.22445414690633522 rtol = 1e-8
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ -4.695916382523908 rtol = 1e-8
    end
end

@testset "CKSAFT" begin
    test_system = system(["carbon dioxide","2-propanol"], "CKSAFT")
    Temp        = 298.15
    Vol         = 1e-4
    z           = create_z(test_system, [0.5,0.5])
    @test typeof(test_system) == OpenSAFT.CKSAFT
    @testset "params" begin
        @test test_system.params.segment[Set(["2-propanol"])] ≈ 3.249 rtol = 1e-8
        @test test_system.params.sigma[Set(["2-propanol"])] ≈ 3.0430915331057595e-10 rtol = 1e-8
        @test test_system.params.sigma[Set(["carbon dioxide", "2-propanol"])] ≈ 3.1070590057871705e-10 rtol = 1e-8
        @test test_system.params.epsilon[Set(["2-propanol"])] ≈ 202.94 rtol = 1e-8
        @test test_system.params.epsilon[Set(["carbon dioxide", "2-propanol"])] ≈ 198.51779876711507 rtol = 1e-8
        @test test_system.params.epsilon_assoc[Set([(Set(["2-propanol"]), "H"), (Set(["2-propanol"]), "e")])] ≈ 2670.0 rtol = 1e-8
        @test test_system.params.bond_vol[Set([(Set(["2-propanol"]), "H"), (Set(["2-propanol"]), "e")])] ≈ 0.021 rtol = 1e-8
        @test test_system.params.n_sites[Set(["2-propanol"])]["e"] ≈ 1 rtol = 1e-8
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_seg(test_system, z, Vol, Temp) ≈ -1.2395529662948277 rtol = 1e-8
        @test OpenSAFT.a_chain(test_system, z, Vol, Temp) ≈ -0.7747586154084931 rtol = 1e-8
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ -1.7937079004096872 rtol = 1e-8
    end
end

@testset "SAFTVRSW" begin
    test_system = system(["water","ethane"], "SAFTVRSW")
    Temp        = 298.15
    Vol         = 1e-4
    z           = create_z(test_system, [0.5, 0.5])
    @test typeof(test_system) == OpenSAFT.SAFTVRSW
    @testset "params" begin
        @test test_system.params.segment[Set(["ethane"])] ≈ 1.33 rtol = 1e-8
        @test test_system.params.sigma[Set(["water"])] ≈ 3.036e-10 rtol = 1e-8
        @test test_system.params.sigma[Set(["water","ethane"])] ≈ 3.6345e-10 rtol = 1e-8
        @test test_system.params.epsilon[Set(["water"])] ≈ 253.3 rtol = 1e-8
        @test test_system.params.epsilon[Set(["water","ethane"])] ≈ 238.6248939234966 rtol = 1e-8
        @test test_system.params.lambda[Set(["water"])] ≈ 1.8 rtol = 1e-8
        @test test_system.params.lambda[Set(["water","ethane"])] ≈ 1.595600082542303 rtol = 1e-8
        @test test_system.params.epsilon_assoc[Set([(Set(["water"]), "H"), (Set(["water"]), "e")])] ≈ 1366.0 rtol = 1e-8
        @test test_system.params.bond_vol[Set([(Set(["water"]), "H"), (Set(["water"]), "e")])] ≈ 1.028e-30 rtol = 1e-8
        @test test_system.params.n_sites[Set(["water"])]["e"] ≈ 2 rtol = 1e-8
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_mono(test_system, z, Vol, Temp) ≈ -1.4659622048306407 rtol = 1e-8
        @test OpenSAFT.a_chain(test_system, z, Vol, Temp) ≈ 0.022703334973543182 rtol = 1e-8
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ -0.5091186885233859 rtol = 1e-8
    end
end

@testset "softSAFT" begin
    test_system = system(["methanol"], "softSAFT")
    Temp        = 298.15
    Vol         = 1e-4
    z           = create_z(test_system, [1.])
    @test typeof(test_system) == OpenSAFT.softSAFT
    @testset "params" begin
        @test test_system.params.segment[Set(["methanol"])] ≈ 1.491 rtol = 1e-8
        @test test_system.params.sigma[Set(["methanol"])] ≈ 3.375e-10 rtol = 1e-8
        @test test_system.params.epsilon[Set(["methanol"])] ≈ 220.4 rtol = 1e-8
        @test test_system.params.epsilon_assoc[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 3213.0 rtol = 1e-8
        @test test_system.params.bond_vol[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 4.847e-27 rtol = 1e-8
        @test test_system.params.n_sites[Set(["methanol"])]["e"] ≈ 1 rtol = 1e-8
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_LJ(test_system, z, Vol, Temp) ≈ -1.299697047509115 rtol = 1e-8
        @test OpenSAFT.a_chain(test_system, z, Vol, Temp) ≈ 0.33553325545339646 rtol = 1e-8
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ -3.7490459421490447 rtol = 1e-8
    end
end

@testset "PCSAFT" begin
    test_system = system(["butane", "ethanol"], "PCSAFT")
    Temp        = 298.15
    Vol         = 1e-4
    z           = create_z(test_system, [0.5, 0.5])
    @test typeof(test_system) == OpenSAFT.PCSAFT
    @testset "params" begin
        @test test_system.params.segment[Set(["ethanol"])] ≈ 2.3827 rtol = 1e-8
        @test test_system.params.sigma[Set(["ethanol"])] ≈ 3.1771e-10 rtol = 1e-8
        @test test_system.params.sigma[Set(["butane", "ethanol"])] ≈ 3.44285e-10 rtol = 1e-8
        @test test_system.params.epsilon[Set(["ethanol"])] ≈ 198.24 rtol = 1e-8
        @test test_system.params.epsilon[Set(["butane", "ethanol"])] ≈ 204.31368602729677 rtol = 1e-8
        @test test_system.params.epsilon_assoc[Set([(Set(["ethanol"]), "H"), (Set(["ethanol"]), "e")])] ≈ 2653.4 rtol = 1e-8
        @test test_system.params.bond_vol[Set([(Set(["ethanol"]), "H"), (Set(["ethanol"]), "e")])] ≈ 0.032384 rtol = 1e-8
        @test test_system.params.n_sites[Set(["ethanol"])]["e"] ≈ 1 rtol = 1e-8
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_hc(test_system, z, Vol, Temp) ≈ 3.114823074155765 rtol = 1e-8
        @test OpenSAFT.a_disp(test_system, z, Vol, Temp) ≈ -6.090736624622517 rtol = 1e-8
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ -2.121606453473655 rtol = 1e-8
    end
end

@testset "sPCSAFT" begin
    test_system = system(["pentane", "methanol"], "sPCSAFT")
    Temp        = 298.15
    Vol         = 1e-4
    z           = create_z(test_system, [0.5, 0.5])
    @test typeof(test_system) == OpenSAFT.sPCSAFT
    @testset "params" begin
        @test test_system.params.segment[Set(["methanol"])] ≈ 2.7921 rtol = 1e-8
        @test test_system.params.sigma[Set(["methanol"])] ≈ 2.651e-10 rtol = 1e-8
        @test test_system.params.sigma[Set(["pentane", "methanol"])] ≈ 3.21195e-10 rtol = 1e-8
        @test test_system.params.epsilon[Set(["methanol"])] ≈ 186.6 rtol = 1e-8
        @test test_system.params.epsilon[Set(["pentane", "methanol"])] ≈ 207.7063311504972 rtol = 1e-8
        @test test_system.params.epsilon_assoc[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 2090.2 rtol = 1e-8
        @test test_system.params.bond_vol[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 0.146 rtol = 1e-8
        @test test_system.params.n_sites[Set(["methanol"])]["e"] ≈ 1 rtol = 1e-8
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_hc(test_system, z, Vol, Temp) ≈ 3.568650770403549 rtol = 1e-8
        @test OpenSAFT.a_disp(test_system, z, Vol, Temp) ≈ -6.994181358803752 rtol = 1e-8
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ -1.7525112985184315 rtol = 1e-8
    end
end

@testset "SAFTVRMie" begin
    test_system = system(["methanol"], "SAFTVRMie")
    Temp        = 298.15
    Vol         = 1e-4
    z           = create_z(test_system, [1.])
    @test typeof(test_system) == OpenSAFT.SAFTVRMie
    @testset "params" begin
        @test test_system.params.segment[Set(["methanol"])] ≈ 1.67034 rtol = 1e-8
        @test test_system.params.sigma[Set(["methanol"])] ≈ 3.2462e-10 rtol = 1e-8
        @test test_system.params.epsilon[Set(["methanol"])] ≈ 307.69 rtol = 1e-8
        @test test_system.params.lambdaA[Set(["methanol"])] ≈ 6.0 rtol = 1e-8
        @test test_system.params.lambdaR[Set(["methanol"])] ≈ 19.235 rtol = 1e-8
        @test test_system.params.epsilon_assoc[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 2062.1 rtol = 1e-8
        @test test_system.params.bond_vol[Set([(Set(["methanol"]), "H"), (Set(["methanol"]), "e")])] ≈ 1.0657e-28 rtol = 1e-8
        @test test_system.params.n_sites[Set(["methanol"])]["e"] ≈ 2 rtol = 1e-8
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_mono(test_system, z, Vol, Temp) ≈ -1.7176380421592838 rtol = 1e-8
        @test OpenSAFT.a_chain(test_system, z, Vol, Temp) ≈ -0.030259270092795967 rtol = 1e-8
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ -3.1565551121889293 rtol = 1e-8
    end
end

@testset "SAFTVRQMie" begin
    test_system = system(["helium"], "SAFTVRQMie")
    Temp        = 298.15
    Vol         = 1e-4
    z           = create_z(test_system, [1.])
    @test typeof(test_system) == OpenSAFT.SAFTVRQMie
    @testset "params" begin
        @test test_system.params.segment[Set(["helium"])] ≈ 1. rtol = 1e-8
        @test test_system.params.sigma[Set(["helium"])] ≈ 2.549e-10 rtol = 1e-8
        @test test_system.params.epsilon[Set(["helium"])] ≈ 10.952 rtol = 1e-8
        @test test_system.params.lambdaA[Set(["helium"])] ≈ 6.0 rtol = 1e-8
        @test test_system.params.lambdaR[Set(["helium"])] ≈ 13.0 rtol = 1e-8
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_mono(test_system, z, Vol, Temp) ≈ 0.12253715358076675 rtol = 1e-8
    end
end

@testset "SAFTgammaMie" begin
    test_system = system(gc("ethanol"), "SAFTgammaMie")
    Temp        = 298.15
    Vol         = exp10(-3.5)
    z           = create_z(test_system, [1.])
    @test typeof(test_system) == OpenSAFT.SAFTgammaMie
    @testset "params" begin
        @test test_system.params.segment[Set(["CH2OH"])] ≈ 2.0 rtol = 1e-8
        @test test_system.params.shapefactor[Set(["CH2OH"])] ≈ 0.58538 rtol = 1e-8
        @test test_system.params.sigma[Set(["CH2OH"])] ≈ 3.4054000000000005e-10 rtol = 1e-8
        @test test_system.params.sigma[Set(["CH2OH","CH3"])] ≈ 3.7413000000000003e-10 rtol = 1e-8
        @test test_system.params.epsilon[Set(["CH2OH"])] ≈ 407.22 rtol = 1e-8
        @test test_system.params.epsilon[Set(["CH2OH","CH3"])] ≈ 333.2 rtol = 1e-8
        @test test_system.params.lambda_a[Set(["CH2OH"])] ≈ 6.0 rtol = 1e-8
        @test test_system.params.lambda_a[Set(["CH2OH","CH3"])] ≈ 6.0 rtol = 1e-8
        @test test_system.params.lambda_r[Set(["CH2OH"])] ≈ 22.699 rtol = 1e-8
        @test test_system.params.lambda_r[Set(["CH2OH","CH3"])] ≈ 18.406912409694556 rtol = 1e-8
        @test test_system.params.epsilon_assoc[Set([(Set(["CH2OH"]), "H"), (Set(["CH2OH"]), "e1")])] ≈ 2097.9 rtol = 1e-8
        @test test_system.params.bond_vol[Set([(Set(["CH2OH"]), "H"), (Set(["CH2OH"]), "e1")])] ≈ 6.2309e-29 rtol = 1e-8
        @test test_system.params.n_sites[Set(["CH2OH"])]["e1"] ≈ 2 rtol = 1e-8
    end
    # We can try to extract specific terms, although I don't know how useful it will be
    @testset "equations" begin
        @test OpenSAFT.a_mono(test_system, z, Vol, Temp) ≈ -1.151043781769667 rtol = 1e-8
        @test OpenSAFT.a_chain(test_system, z, Vol, Temp) ≈ -0.1255227354789658 rtol = 1e-8
        @test OpenSAFT.a_assoc(test_system, z, Vol, Temp) ≈ -1.9386416653191778 rtol = 1e-8
    end
end
