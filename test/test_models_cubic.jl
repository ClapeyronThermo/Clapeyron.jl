@testset "Cubic models" begin
    @printline
    let T = 333.15, V = 1e-3,p = 1e5,z = [0.5,0.5],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "vdW Models" begin
        @testset "Default vdW" begin
            system = vdW(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.7088380780265725 rtol = 1e-6
            @test Clapeyron.cubic_poly(system, p, T, z)[1][1] ≈ -0.0002475728429728521 rtol = 1e-6
            @test Clapeyron.cubic_p(system, V, T, z) ≈ Clapeyron.pressure(system, V, T, z) rtol = 1e-6
        end

        @testset "Clausius" begin
            system = Clausius(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2945136000972637 rtol = 1e-6
        end

        @testset "Berthelot" begin
            system = Berthelot(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.5282310106210891 rtol = 1e-6
        end
    end

    @testset "RK Models" begin
        @testset "Default RK" begin
            system = RK(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.9825375012134132 rtol = 1e-6
            @test Clapeyron.cubic_poly(system, p, T, z)[1][1] ≈ -0.00022230043592123767 rtol = 1e-6
            @test Clapeyron.cubic_p(system, V, T, z) ≈ Clapeyron.pressure(system, V, T, z) rtol = 1e-6
        end

        @testset "SRK" begin
            system = SRK(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2572506872856557 rtol = 1e-6
        end

        @testset "PSRK" begin
            system = PSRK(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2265133881057408 rtol = 1e-6
        end

        @testset "tcRK" begin
            system = tcRK(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2519495469149748 rtol = 1e-6
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
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.5729126903890258 rtol = 1e-6
        end
    end

    @testset "PR Models" begin
        @testset "Default PR" begin
            system = PR(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.244774062489359 rtol = 1e-6
            @test Clapeyron.cubic_poly(system, p, T, z)[1][1] ≈ -0.0002328543992909459 rtol = 1e-6
            @test Clapeyron.cubic_p(system, V, T, z) ≈ Clapeyron.pressure(system, V, T, z) rtol = 1e-6
        end

        @testset "PR78" begin
            system = PR78(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.246271020258271 rtol = 1e-6
        end

        @testset "VTPR" begin
            system = VTPR(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.234012667532541 rtol = 1e-6
        end

        @testset "UMRPR" begin
            system = UMRPR(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.1447330557895619 rtol = 1e-6
        end

        @testset "QCPR" begin
            system = QCPR(["neon","helium"])
            @test Clapeyron.a_res(system, V, 25, z) ≈ -0.04727884027682022 rtol = 1e-6
            @test Clapeyron.lb_volume(system,25,z) ≈ 1.3601716423130568e-5 rtol = 1e-6
            _a,_b,_c = Clapeyron.cubic_ab(system,V,25,z)
            @test _a ≈ 0.012772722389495079 rtol = 1e-6
            @test _b ≈ 1.0728356231510917e-5 rtol = 1e-6
            @test _c ≈ -2.87335e-6 rtol = 1e-6
            #test for the single component branch
            system1 = QCPR(["helium"])
            a1 = Clapeyron.a_res(system, V, 25, [0.0,1.0])
            @test Clapeyron.a_res(system1, V, 25, [1.0])  ≈ a1
        end

        @testset "cPR" begin
            system = cPR(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2438144556398565 rtol = 1e-6
        end

        @testset "tcPR" begin
            system = tcPR(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.254190142912733 rtol = 1e-6
        end

        @testset "tcPR + Wilson (Res)" begin
            system = tcPRW(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2106271631685903 rtol = 1e-6
        end

        @testset "EPPR78" begin
           system = EPPR78(["benzene","isooctane"])
           @test Clapeyron.a_res(system, V, T, z) ≈ -1.1415931771485368 rtol = 1e-6
        end

        @testset "PR w/ BMAlpha" begin
            system = PR(["ethane","undecane"];alpha = BMAlpha)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2445088818575114 rtol = 1e-6
        end

        @testset "PR w/ TwuAlpha" begin
            system = PR(["ethane","undecane"];alpha = TwuAlpha)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2650756692036234 rtol = 1e-6
        end

        @testset "PR w/ MTAlpha" begin
            system = PR(["ethane","undecane"];alpha = MTAlpha)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2542346631395425 rtol = 1e-6
        end

        @testset "PR w/ LeiboviciAlpha" begin
            system = PR(["ethane","undecane"];alpha = LeiboviciAlpha)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2480909069722526 rtol = 1e-6
        end

        @testset "PR w/ RackettTranslation" begin
            system = PR(["ethane","undecane"];translation = RackettTranslation)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2453853855058576 rtol = 1e-6
        end

        @testset "PR w/ MTTranslation" begin
            system = PR(["ethane","undecane"];translation = MTTranslation)
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.243918482158021 rtol = 1e-6
        end

        @testset "PR w/ HVRule" begin
            system = PR(["methanol","benzene"];mixing = HVRule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.6329827794751909 rtol = 1e-6
        end

        @testset "PR w/ MHV1Rule" begin
            system = PR(["methanol","benzene"];mixing = MHV1Rule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.6211441939544694 rtol = 1e-6
        end

        @testset "PR w/ MHV2Rule" begin
            system = PR(["methanol","benzene"];mixing = MHV2Rule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.6210843663212396 rtol = 1e-6
        end

        @testset "PR w/ LCVMRule" begin
            system = PR(["methanol","benzene"];mixing = LCVMRule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.6286241404575419 rtol = 1e-6
        end

        @testset "PR w/ WSRule" begin
            system = PR(["methanol","benzene"];mixing = WSRule, activity=Wilson)
            @test Clapeyron.a_res(system, V, T, z) ≈ -0.6690864227574802 rtol = 1e-6
        end
    end

    @testset "KU Models" begin
        system = KU(["ethane","undecane"])
        @test Clapeyron.a_res(system, V, T, z) ≈ -1.2261554720898895 rtol = 1e-6
        @test Clapeyron.cubic_p(system, V, T, z) ≈ Clapeyron.pressure(system, V, T, z) rtol = 1e-6
    end

    @testset "RKPR Models" begin
        system = RKPR(["ethane","undecane"])
        @test Clapeyron.a_res(system, V, T, z) ≈ -1.2877492838069213 rtol = 1e-6
    end

    @testset "Patel Teja Models" begin
        @testset "Patel Teja" begin
            system = PatelTeja(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2284322450064429 rtol = 1e-6
            @test Clapeyron.cubic_p(system, V, T, z) ≈ Clapeyron.pressure(system, V, T, z) rtol = 1e-6
        end
        @testset "PTV" begin
            system = PTV(["ethane","undecane"])
            @test Clapeyron.a_res(system, V, T, z) ≈ -1.2696422558756286 rtol = 1e-6
            @test Clapeyron.cubic_p(system, V, T, z) ≈ Clapeyron.pressure(system, V, T, z) rtol = 1e-6
        end
    end
    end
    @printline
end
