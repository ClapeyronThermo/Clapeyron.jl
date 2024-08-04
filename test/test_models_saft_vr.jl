GC.gc()

@testset "SAFT-VR-Mie Models" begin
    @printline
    let T = 298.15, V = 1e-4,z1 = Clapeyron.SA[1.0],z = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "SAFTVRMie" begin
        system = SAFTVRMie(["methanol", "water"])
        @test Clapeyron.a_mono(system, V, T, z) ≈ -0.9729134860869052 rtol = 1e-6
        _a_chain = Clapeyron.a_chain(system, V, T, z)
        _a_disp  = Clapeyron.a_disp(system, V, T, z)
        @test _a_chain ≈ -0.02834738013535014 rtol = 1e-6
        @test Clapeyron.a_dispchain(system,V,T,z) - _a_chain ≈ _a_disp rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -4.180807072390184 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
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

    @testset "SAFTgammaMie" begin
        system = SAFTgammaMie(["methanol","butane"])
        V_γMie = exp10(-3.5)
        @test Clapeyron.a_mono(system, V_γMie, T, z) ≈ -1.0400249396482548 rtol = 1e-6
        @test Clapeyron.a_chain(system, V_γMie, T, z) ≈ -0.07550931466871749 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V_γMie, T, z) ≈ -0.8205840455850311 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end
    
    @testset "SAFTVRMieGV" begin # SSTODO: test these
        system1 = SAFTVRMieGV(["carbon dioxide", "acetone"])
        # system2 = SAFTVRMieGV(["carbon dioxide", "water"])
        let T = 298.15, V = 4.9914e-5, z = [0.5, 0.5]
        @test Clapeyron.a_mp(system1, V, T, z) ≈ -2.25103330151645    rtol = 1e-6 
        # @test Clapeyron.a_mp(system2, V, T, z) ≈ -0.1392358363758833 rtol = 1e-6
        test_gibbs_duhem(system1,V,T,z)
        # test_gibbs_duhem(system2,V,T,z)
        GC.gc()
        end
    end
    
    @testset "structSAFTgammaMie" begin
        species = [("ethanol",["CH3"=>1,"CH2OH"=>1],[("CH3","CH2OH")=>1]),
                   ("octane",["CH3"=>2,"CH2"=>6],[("CH3","CH2")=>2,("CH2","CH2")=>5])]

        system = structSAFTgammaMie(species)
        V_γMie = exp10(-3.5)
        @test Clapeyron.a_chain(system, V_γMie, T, z) ≈ -0.11160851237651681 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z,rtol = 1e-12)
        GC.gc()
    end
    end
    @printline
end
