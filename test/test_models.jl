using Clapeyron, Test

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

    @testset "LJSAFT" begin
        system = LJSAFT(["ethane","1-propanol"])
        z = [0.5, 0.5]
        @test Clapeyron.a_seg(system, V, T, z) ≈ -2.207632433058473 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ -0.04577483379871112 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.3009761155167205 rtol = 1e-6
    end

    @testset "CPA" begin
        system = CPA(["ethanol","benzene"])
        z = [0.5, 0.5]
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.1575210505284332 rtol = 1e-6
    end

    @testset "sCPA" begin
        system = sCPA(["water","carbon dioxide"])
        z = [0.5, 0.5]
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.957518287413705 rtol = 1e-6
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
    end

    @testset "SAFTVRMie" begin
        system = SAFTVRMie(["methanol", "water"])
        z = [0.5, 0.5]
        @test Clapeyron.a_mono(system, V, T, z) ≈ -0.9729134860869052 rtol = 1e-6
        _a_chain = Clapeyron.a_chain(system, V, T, z)
        _a_disp  =  Clapeyron.a_disp(system, V, T, z)
        @test _a_chain ≈ -0.02834738013535014 rtol = 1e-6
        @test Clapeyron.a_dispchain(system,V,T,z) - _a_chain ≈ _a_disp rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -4.180807072390184 rtol = 1e-6
    end

    @testset "SAFTVRQMie" begin
        system = SAFTVRQMie(["helium"])
        z = [1.]
        @test Clapeyron.a_mono(system, V, T, z) ≈ 0.12286776703976324 rtol = 1e-6
    end

    @testset "SAFTgammaMie" begin
        system = SAFTgammaMie(["methanol","butane"])
        V_γMie = exp10(-3.5)
        z = [0.5,0.5]
        @test Clapeyron.a_mono(system, V_γMie, T, z) ≈ -1.0400249396482548 rtol = 1e-6
        @test Clapeyron.a_chain(system, V_γMie, T, z) ≈ -0.07550931466871749 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V_γMie, T, z) ≈ -0.8205840455850311 rtol = 1e-6
    end
end

@testset "Cubic models" begin
    T = 333.15
    V = 1e-3
    p = 1e5
    z = [0.5, 0.5]

    @testset "vdW Models" begin
            system = vdW(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.7088380780265725 rtol = 1e-6
            @test Clapeyron.cubic_poly(system, p, T, z)[1][1] ≈ -0.0002475728429728521 rtol = 1e-6
            @test Clapeyron.cubic_abp(system, V, T, z)[1] ≈ 2.591962221523441 rtol = 1e-6
    end

    @testset "RK Models" begin
        @testset "Default RK" begin
            system = RK(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.9825375012134132 rtol = 1e-6
            @test Clapeyron.cubic_poly(system, p, T, z)[1][1] ≈ -0.00022230043592123767 rtol = 1e-6
            @test Clapeyron.cubic_abp(system, V, T, z)[1] ≈ 3.357807485319651 rtol = 1e-6
        end

        @testset "SRK" begin
            system = SRK(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2572506872856557 rtol = 1e-6
        end

        @testset "PSRK" begin
            system = PSRK(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2265133881057408 rtol = 1e-6
        end

        @testset "RK w/ BMAlpha" begin
            system = RK(["ethane","undecane"];alpha = BMAlpha)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2569334957019538 rtol = 1e-6
        end

        @testset "RK w/ PenelouxTranslation" begin
            system = RK(["ethane","undecane"];translation = PenelouxTranslation)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.9800472116681871 rtol = 1e-6
        end

        @testset "RK w/ KayRule" begin
            system = RK(["ethane","undecane"];mixing = KayRule)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.8176850121211936 rtol = 1e-6
        end

        @testset "RK w/ HVRule" begin
            system = RK(["methanol","benzene"];mixing = HVRule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.5209112693371991 rtol = 1e-6
        end

        @testset "RK w/ MHV1Rule" begin
            system = RK(["methanol","benzene"];mixing = MHV1Rule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.5091065987876959 rtol = 1e-6
        end

        @testset "RK w/ MHV2Rule" begin
            system = RK(["methanol","benzene"];mixing = MHV2Rule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.5071048453490448 rtol = 1e-6
        end

        @testset "RK w/ WSRule" begin
            system = RK(["methanol","benzene"];mixing = WSRule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.572912690389026 rtol = 1e-6
        end
    end

    @testset "PR Models" begin
        @testset "Default PR" begin
            system = PR(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.244772730766631 rtol = 1e-6
            @test Clapeyron.cubic_poly(system, p, T, z)[1][1] ≈ -0.00023285390449318037 rtol = 1e-6
            @test Clapeyron.cubic_abp(system, V, T, z)[1] ≈ 4.268630968024985 rtol = 1e-6
        end

        @testset "PR78" begin
            system = PR78(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.246269686941749 rtol = 1e-6
        end

        @testset "VTPR" begin
            system = VTPR(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.23401133834732 rtol = 1e-6
        end

        @testset "UMRPR" begin
            system = UMRPR(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.1447318247939071 rtol = 1e-6
        end

        @testset "PR w/ BMAlpha" begin
            system = PR(["ethane","undecane"];alpha = BMAlpha)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.244507550417118 rtol = 1e-6
        end

        @testset "PR w/ TwuAlpha" begin
            system = PR(["ethane","undecane"];alpha = TwuAlpha)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2650743158660063 rtol = 1e-6
        end

        @testset "PR w/ MTAlpha" begin
            system = PR(["ethane","undecane"];alpha = MTAlpha)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2542333213442207 rtol = 1e-6
        end
        
        @testset "PR w/ RackettTranslation" begin
            system = PR(["ethane","undecane"];translation = RackettTranslation)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2453840474643165 rtol = 1e-6
        end

        @testset "PR w/ MTTranslation" begin
            system = PR(["ethane","undecane"];translation = MTTranslation)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.24391715980424 rtol = 1e-6
        end

        @testset "PR w/ HVRule" begin
            system = PR(["methanol","benzene"];mixing = HVRule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.632982061564318 rtol = 1e-6
        end

        @testset "PR w/ MHV1Rule" begin
            system = PR(["methanol","benzene"];mixing = MHV1Rule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.6211434867194446 rtol = 1e-6
        end

        @testset "PR w/ MHV2Rule" begin
            system = PR(["methanol","benzene"];mixing = MHV2Rule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.6210836570941939 rtol = 1e-6
        end

        @testset "PR w/ LCVMRule" begin
            system = PR(["methanol","benzene"];mixing = LCVMRule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.6286234264772204 rtol = 1e-6
        end

        @testset "PR w/ WSRule" begin
            system = PR(["methanol","benzene"];mixing = WSRule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.669085674824878 rtol = 1e-6
        end
    end
end

@testset "Activity models" begin
    T = 333.15
    p = 1e5
    z = [0.5,0.5]

    @testset "Wilson" begin
        system = Wilson(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.530046633499114 rtol = 1e-6
    end

    @testset "NRTL" begin
        system = NRTL(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5309354738922405 rtol = 1e-6
    end

    @testset "UNIQUAC" begin
        system = UNIQUAC(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.3630421218486388 rtol = 1e-6
    end

    @testset "UNIFAC" begin
        system = UNIFAC(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5322232657797463 rtol = 1e-6
    end

    @testset "ogUNIFAC" begin
        system = ogUNIFAC(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5133696314734384 rtol = 1e-6
    end

    @testset "COSMOSAC02" begin
        system = COSMOSAC02(["water","ethanol"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.3871817962565904 rtol = 1e-6
    end

    @testset "COSMOSAC10" begin
        system = COSMOSAC10(["water","ethanol"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.4015660588643404 rtol = 1e-6
    end

    @testset "COSMOSACdsp" begin
        system = COSMOSACdsp(["water","ethanol"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.4398951117248127 rtol = 1e-6
    end
end

@testset "Ideal models" begin
    T = 298.15
    V = 1e-4
    z = [1.]

    @testset "Joback" begin
        system = JobackIdeal(["hexane"])
        @test Clapeyron.C_p(system,298.15) ≈ 143.22076150138616 rtol = 1e-6
        @test Clapeyron.T_b(system) ≈ 336.88 rtol = 1e-6
        @test Clapeyron.crit_pure(system)[1] ≈ 500.2728274871347 rtol = 1e-6
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 9.210841420941021 rtol = 1e-6
    end

    @testset "Reid" begin
        system = ReidIdeal(["butane"])
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 9.210842104089576 rtol = 1e-6
    end

    @testset "Walker" begin
        system = WalkerIdeal(["hexane"])
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 179.51502015696653 rtol = 1e-6
    end

    @testset "Monomer" begin
        system = MonomerIdeal(["hexane"])
        @test Clapeyron.a_ideal(system,V,T,z) ≈ -20.368750666236373 rtol = 1e-6
    end
end

@testset "Multi-parameter models" begin
    T = 298.15
    V = 1e-4

    @testset "IAPWS95" begin
        z = [1.]
        system = IAPWS95()
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 7.932118505056652 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -2.1152657050144347e14 rtol = 1e-6
    end

    @testset "PropaneRef" begin
        z   = [1.]
        system = PropaneRef()
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 0.6426994942361217 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -2.436280448227229 rtol = 1e-6
    end

    @testset "GERG2008" begin
        z   = [1.]
        system = GERG2008(["water"])
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 4.500099760879548 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -10.122119572808764 rtol = 1e-6
        z   = [0.25,0.25,0.25,0.25]
        system = GERG2008(["water","carbon dioxide","hydrogen sulfide","argon"])
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 3.1135835641766594 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -1.1706377677539772 rtol = 1e-6
    end
end

@testset "SPUNG models" begin
    T = 298.15
    V = 1e-4
    z = [1.]

    @testset "SRK" begin
        system = SPUNG(["ethane"])
        @test Clapeyron.shape_factors(system, V, T, z)[1] ≈ 0.7289071312821193 rtol = 1e-6
    end

    @testset "PCSAFT" begin
        system = SPUNG(["ethane"],PropaneRef(),PCSAFT(["ethane"]),PCSAFT(["propane"]))
        @test Clapeyron.shape_factors(system, V, T, z)[1] ≈ 1.3499576779924594 rtol = 1e-6
    end
end

@testset "lattice models" begin
    T = 298.15
    V = 1e-4
    z = [1.]
    system = Clapeyron.SanchezLacombe(["carbon dioxide"])
    @test Clapeyron.a_res(system, V, T, z) ≈ -0.9511044462267396 rtol = 1e-6
end
