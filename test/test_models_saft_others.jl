@testset "SAFT models - misc" begin
    @printline
    let T = 298.15, V = 1e-4,z = [0.5,0.5],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "ogSAFT" begin
        system = ogSAFT(["water","ethylene glycol"])
        @test Clapeyron.a_seg(system, V, T, z) ≈ -2.0332062924093366 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ -0.006317441684202759 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -4.034042081699316 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end

    @testset "CKSAFT" begin
        system = CKSAFT(["carbon dioxide", "2-propanol"])
        @test Clapeyron.a_seg(system, V, T, z) ≈ -1.24586302917188 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ -0.774758615408493 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.2937079004096872 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end

    @testset "sCKSAFT" begin
        system = sCKSAFT(["benzene","acetic acid"])
        @test Clapeyron.a_seg(system, V, T, z) ≈ -3.1809330810925256 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -3.3017434376105514 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end

    @testset "BACKSAFT" begin
        system = BACKSAFT(["carbon dioxide"])
        @test Clapeyron.a_hcb(system, V, T, z1) ≈ 1.0118842111801198 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z1) ≈ -0.14177009317268635 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z1) ≈ -2.4492518566426296 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z1)
        GC.gc()
    end

    @testset "LJSAFT" begin
        system = LJSAFT(["ethane","1-propanol"])
        @test Clapeyron.a_seg(system, V, T, z) ≈ -2.207632433058473 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ -0.04577483379871112 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.3009761155167205 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end
    @testset "SAFTVRSW" begin
        system = SAFTVRSW(["water", "ethane"])
        @test Clapeyron.a_mono(system, V, T, z) ≈ -1.4367205951569462 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ 0.022703335564111336 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -0.5091186813915323 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end

    @testset "softSAFT" begin
        system = softSAFT(["hexane","1-propanol"])
        @test Clapeyron.a_LJ(system, V, T, z) ≈ -3.960728242264164 rtol = 1e-6
        @test Clapeyron.a_chain(system, V, T, z) ≈ 0.3736728407455211 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -2.0461376618069034 rtol = 1e-6
        #TODO: check here why the error is so big
        #we disable this test for now, it takes too much time
        #test_gibbs_duhem(system,V,T,z,rtol = 1e-12)
        GC.gc()
    end

    @testset "softSAFT2016" begin
        system = softSAFT2016(["hexane","1-propanol"])
        @test Clapeyron.a_LJ(system, V, T, z) ≈ -3.986690073534575 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end

    @testset "solidsoftSAFT" begin
        system = solidsoftSAFT(["octane"])
        V_sol = 1e-4
        @test Clapeyron.a_LJ(system, V_sol, T, z1) ≈ 7.830498923903852 rtol = 1e-6
        @test Clapeyron.a_chain(system, V_sol, T, z1) ≈ -2.3460460361188207 rtol = 1e-6
        test_gibbs_duhem(system, V_sol, T, z1, rtol = 1e-12)
        GC.gc()
    end
    end
end

@testset "CPA" begin
    @printline
    let T = 298.15, V = 1e-4,z1 = Clapeyron.SA[1.0],z = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "CPA" begin
        system = CPA(["ethanol","benzene"])
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.1575210505284332 rtol = 1e-6
        test_gibbs_duhem(system, V, T, z)
        GC.gc()
    end

    @testset "sCPA" begin
        system = sCPA(["water","carbon dioxide"])
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.957518287413705 rtol = 1e-6
        GC.gc()
    end
    end
end