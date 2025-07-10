@testset "Electrolyte models" begin
    @printline

    let T = 298.15, V = 1e-3, m = 1, w = [0.99,0.005,0.005],Z = [0,1,-1];

    @testset "ConstRSP" begin
        system = ConstRSP()
        @test Clapeyron.dielectric_constant(system, V, T, w, Z) ≈ 78.38484961 rtol = 1e-6
    end

    @testset "ZuoFurst" begin
        system = ZuoFurst()
        @test Clapeyron.dielectric_constant(system, V, T, w, Z) ≈ 78.30270731614397 rtol = 1e-6
    end

    @testset "LinMixRSP" begin
        system = LinMixRSP(["water"],["sodium","chloride"])
        @test Clapeyron.dielectric_constant(system, V, T, w, Z) ≈ 77.68100111390001 rtol = 1e-6
    end

    @testset "Schreckenberg" begin
        system = Schreckenberg(["water"],["sodium","chloride"])
        @test Clapeyron.dielectric_constant(system, 1.8e-5 , T, w, Z) ≈ 77.9800485493879 rtol = 1e-6
    end 

    @testset "ePCSAFT" begin
        system = ePCSAFT(["water"],["sodium","chloride"])
        salts = [("sodium chloride",["sodium"=>1,"chloride"=>1])]
        z = molality_to_composition(system,salts,m)
        #this result is using normal water parameters. From Clapeyron 0.6.9 onwards, "water" will use T-dependend sigma on pharmaPCSAFT.
        #@test Clapeyron.a_res(system, V, T, z) ≈ -1.1094732499161208 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -1.0037641675166717 rtol = 1e-6
    end

    @testset "eSAFTVRMie" begin
        system = eSAFTVRMie(["water"],["sodium","chloride"])
        salts = [("sodium chloride",["sodium"=>1,"chloride"=>1])]
        z = molality_to_composition(system,salts,m)
        @test Clapeyron.a_res(system, V, T, z) ≈ -6.8066182496746634 rtol = 1e-6
    end

    @testset "SAFTVREMie" begin
        system = SAFTVREMie(["water"],["sodium","chloride"])
        salts = [("sodium chloride",["sodium"=>1,"chloride"=>1])]
        z = molality_to_composition(system,salts,m)
        @test Clapeyron.a_res(system, V, T, z) ≈ -5.37837866742139 rtol = 1e-6
    end

    @testset "SAFTVREMie + MSA" begin
        system = SAFTVREMie(["water"],["sodium","chloride"], ionmodel = MSA)
        salts = [("sodium chloride",["sodium"=>1,"chloride"=>1])]
        z = molality_to_composition(system,salts,m)
        @test Clapeyron.a_res(system, V, T, z) ≈ -2.266066359288185 rtol = 1e-6
    end

    @testset "SAFTVREMie + Born" begin
        system = SAFTVREMie(["water"],["sodium","chloride"], ionmodel = Born)
        salts = [("sodium chloride",["sodium"=>1,"chloride"=>1])]
        z = molality_to_composition(system,salts,m)
        @test Clapeyron.a_res(system, V, T, z) ≈ -4.907123437550198 rtol = 1e-6
    end

    @testset "SAFTVREMie + MSAID" begin
        ionmodel = MSAID(["water"],["sodium","chloride"];
                    userlocations = (;
                    charge = [0,1,-1],
                    sigma  = [4.,4.,4.],
                    dipole = [2.2203,0.,0.]))

        system = SAFTVREMie(["water"],["sodium","chloride"], ionmodel = ionmodel)
        iondata = Clapeyron.iondata(system,1.8e-5,298.15, [0.9,0.05,0.05])
        @test Clapeyron.a_res(ionmodel, 1.8e-5,298.15, [0.9,0.05,0.05],iondata) ≈ -224.13923847423726 rtol = 1e-6
    end

    @testset "SAFTgammaEMie" begin
        system = SAFTgammaEMie(["water"],["sodium","acetate"])
        test_scales(system)
        test_recombine(system)
        salts = [("sodium acetate",["sodium"=>1,"acetate"=>1])]
        z = molality_to_composition(system,salts,m)
        #@test Clapeyron.a_res(system, V, T, z) ≈ -5.090195753415607 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -5.401187603151029 rtol = 1e-6
        
        model = SAFTgammaEMie(["water"],["calcium","chloride"]) #issue 349
        @test model isa EoSModel
    end
    @printline
end
end
