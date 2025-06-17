GC.gc()

@testset "SAFT-VR-Mie Models" begin
    @printline
    let T = 298.15, V = 1e-4,z1 = Clapeyron.SA[1.0],z = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "SAFTVRMie" begin
        system = SAFTVRMie(["methanol", "water"])
        @test Clapeyron.a_mono(system, V, T, z) ≈ -0.9729139704318698 rtol = 1e-6
        _a_chain = Clapeyron.a_chain(system, V, T, z)
        _a_disp  = Clapeyron.a_disp(system, V, T, z)
        @test _a_chain ≈ -0.028347378889242814 rtol = 1e-6
        @test Clapeyron.a_dispchain(system,V,T,z) - _a_chain ≈ _a_disp rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -4.18080707238976 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        test_recombine(system)
        GC.gc()
    end

    @testset "SAFTVRMieGV" begin
        system = SAFTVRMieGV(["benzene","acetone"])
        V_GV = 8e-5
        @test Clapeyron.a_mp(system, V_GV, T, z) ≈ -0.7521858819355216 rtol = 1e-6
        test_gibbs_duhem(system,V_GV,T,z)
        GC.gc()
        test_recombine(system)
    end

    @testset "SAFTVRQMie" begin
        system = SAFTVRQMie(["helium"])
        @test Clapeyron.a_mono(system, V, T, z1) ≈ 0.12286776703976324 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z1)
        GC.gc()
    end

    @testset "SAFTVRSMie" begin
        system = SAFTVRSMie(["carbon dioxide"])
        V_sol = 3e-5
        @test Clapeyron.a_mono(system, V_sol, T, z1) ≈ 0.43643302846919896 rtol = 1e-6
        @test Clapeyron.a_chain(system, V_sol, T, z1) ≈ -0.4261294644079463 rtol = 1e-6
        test_gibbs_duhem(system,V_sol,T,z1,rtol = 1e-12)
        GC.gc()
    end

    @testset "SAFTVRMie15" begin
        v15 = 1/(1000*1000*1.011/18.015)
        T15 = 290.0
        vr15 = SAFTVRMie15("water")
        #Dufal, table 4, 290K, f_OH(free) = 0.089
        @test Clapeyron.X(vr15,v15,T15,z1)[1][1] ≈ 0.08922902098124778 rtol = 1e-6
    end

    @testset "SAFTgammaMie" begin
        system = SAFTgammaMie(["methanol","butane"])
        V_γMie = exp10(-3.5)
        @test Clapeyron.a_mono(system, V_γMie, T, z) ≈ -1.0400249396482548 rtol = 1e-6
        @test Clapeyron.a_chain(system, V_γMie, T, z) ≈ -0.07550931466871749 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V_γMie, T, z) ≈ -0.8205840455850311 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        test_scales(system)
        test_recombine(system)
        GC.gc()
    end
    
    @testset "structSAFTgammaMie" begin
        species = [("ethanol",["CH3"=>1,"CH2OH"=>1],[("CH3","CH2OH")=>1]),
                   ("octane",["CH3"=>2,"CH2"=>6],[("CH3","CH2")=>2,("CH2","CH2")=>5])]

        system = structSAFTgammaMie(species)
        V_γMie = exp10(-3.5)
        @test Clapeyron.a_chain(system, V_γMie, T, z) ≈ -0.11160851237651681 rtol = 1e-6
        test_gibbs_duhem(system,V_γMie,T,z,rtol = 1e-12)
        GC.gc()
    end
    end
    @printline
end
