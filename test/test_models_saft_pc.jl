@testset "PCSAFT Models" begin
    @printline
    let T = 298.15, V = 1e-4,z1 = Clapeyron.SA[1.0],z = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "PCSAFT" begin
        system = PCSAFT(["butane", "ethanol"])
        @test Clapeyron.a_hc(system, V, T, z) ≈ 3.1148229872928654 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z) ≈ -6.090736508783152 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.6216064387201956 rtol = 1e-6
        test_gibbs_duhem(system, V, T, z)
        GC.gc()
    end
    @printline
    @testset "PCPSAFT" begin
        system = PCPSAFT(["acetone", "butane", "DMSO"])
        set_k!(system,zeros(3,3))
        set_l!(system,zeros(3,3))
        @test Clapeyron.a_polar(system, V, T, z3) ≈ -0.6541688650413224 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z3)
        GC.gc()
    end
    @printline
    @testset "QPCPSAFT" begin
        system1 = QPCPSAFT(["carbon dioxide", "acetone", "hydrogen sulfide"])
        system2 = QPCPSAFT(["carbon dioxide", "chlorine", "carbon disulfide"])
        @test Clapeyron.a_mp(system1, V, T, z3) ≈ -0.37364363283985724 rtol = 1e-6
        @test Clapeyron.a_mp(system2, V, T, z3) ≈ -0.1392358363758833 rtol = 1e-6
        test_gibbs_duhem(system1,V,T,z3)
        test_gibbs_duhem(system2,V,T,z3)
        GC.gc()
    end
    @printline
    @testset "gcPCPSAFT" begin
        system = gcPCPSAFT(["acetone", "ethane","ethanol"],mixing = :homo)
        test_gibbs_duhem(system,V,T,z3)
        GC.gc()
    end
    @printline
    @testset "sPCSAFT" begin
        system = sPCSAFT(["propane", "methanol"])
        @test Clapeyron.a_hc(system, V, T, z) ≈ 2.024250583187793 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z) ≈ -4.138653131750594 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.1459701721909195 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end
    @printline
    @testset "gcsPCSAFT" begin
        system = gcsPCSAFT(["acetone", "ethane"])
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end
    @printline
    @testset "CP-PCSAFT" begin
        system = CPPCSAFT(["butane", "propane"])
        @test Clapeyron.a_hc(system, V, T, z) ≈ 3.856483933013827 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z) ≈ -6.613302753897683 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end
    @printline
    @testset "GEPCSAFT" begin
        system = GEPCSAFT(["propane", "methanol"])
        @test Clapeyron.a_hc(system, V, T, z) ≈ 1.6473483928460233 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z) ≈ -3.271039575934372 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -1.9511233680313027 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end
    @printline
    @testset "gcPCSAFT" begin
        species = [("ethanol",["CH3"=>1,"CH2"=>1,"OH"=>1],[("CH3","CH2")=>1,("OH","CH2")=>1]),
                   ("hexane",["CH3"=>2,"CH2"=>4],[("CH3","CH2")=>2,("CH2","CH2")=>3])]

        system = gcPCSAFT(species)
        @test Clapeyron.a_hc(system, V, T, z) ≈ 5.485662509904188 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z) ≈ -10.594659479487497 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -0.9528180944200482 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        GC.gc()
    end
    @printline
    @testset "ADPCSAFT" begin
        system = ADPCSAFT(["water"])
        @test Clapeyron.a_hs(system, V, T, z1) ≈ 0.3130578789492178 rtol = 1e-6
        @test Clapeyron.a_disp(system, V, T, z1) ≈ -1.2530666693292463 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z1) ≈ -3.805796041192079 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z1)
        GC.gc()
    end
    end
    @printline
    
    @testset "DAPT" begin
        #try to run dapt in it's own scope
        let T = 298.15, V = 1e-4,z1 = Clapeyron.SA[1.0];
            system = DAPT(["water"])
            @test Clapeyron.a_hs(system, V, T, z1) ≈ 0.35240995905438116 rtol = 1e-6
            @test Clapeyron.a_disp(system, V, T, z1) ≈ -1.7007754776344663 rtol = 1e-6
            @test Clapeyron.a_assoc(system, V, T, z1) ≈ -1.815041612389342 rtol = 1e-6
            test_gibbs_duhem(system,V,T,z1)
            GC.gc()
        end
    end
end
@printline