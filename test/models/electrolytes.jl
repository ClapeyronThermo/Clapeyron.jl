@testset verbose = true "Electrolyte models" begin
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
        test_repr(system,str = ["Neutral Model:","Ion Model:","RSP Model"])
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
        #salts = [("sodium chloride",["sodium"=>1,"chloride"=>1])]
        #z = molality_to_composition(system,salts,m)
        z = molality_to_composition(system,m)
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
        @test Clapeyron.a_res(ionmodel, 1.8e-5,298.15, [0.9,0.05,0.05],iondata) ≈ -224.26174931795842 rtol = 1e-6
    end

    @testset "SAFTgammaEMie" begin
        system = SAFTgammaEMie(["water"],["sodium","acetate"])
        test_scales(system)
        test_recombine(system)
        #salts = [("sodium acetate",["sodium"=>1,"acetate"=>1])]
        #z = molality_to_composition(system,salts,m)
        z = molality_to_composition(system,m)
        #@test Clapeyron.a_res(system, V, T, z) ≈ -5.090195753415607 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -5.401187603151029 rtol = 1e-6
        
        model = SAFTgammaEMie(["water"],["calcium","chloride"]) #issue 349
        @test model isa EoSModel
    end

    @testset "Liquid State Theory" begin
        model = LS([("polymer", ["+"=>10, "-"=>10], [("+","-")=>1, ("+","+")=>4, ("-","-")=>4]),
                   ("counterion", ["-"=>1], [])])
        @test Clapeyron.a_res(model,0.03,300.0,[0.5,0.5]) ≈ -1.5683657474962023 rtol = 1e-6
    end

    @testset "Pitzer" begin
        model = Pitzer(["water"], ["sodium", "potassium", "chloride"]; require_mixing = true)
        test_repr(model, str = ["Pitzer Electrolyte Activity Model", "Temperature Model: constant"])
        @test model isa Clapeyron.ElectrolyteModel && model isa PitzerModel
        @test isempty(pitzer_missing_binary_parameters(model))
        @test isempty(pitzer_missing_mixing_parameters(model).theta)
        @test isempty(pitzer_missing_mixing_parameters(model).psi)

        salts = [("NaCl", ["sodium" => 1, "chloride" => 1]),
                 ("KCl",  ["potassium" => 1, "chloride" => 1])]
        z = molality_to_composition(model, salts, [1.0, 1.0])
        @test composition_to_molality(model, z) ≈ [0.0, 1.0, 1.0, 2.0] rtol = 1e-12
        @test ionic_strength(model, salts, [1.0, 1.0]) ≈ 2.0 rtol = 1e-12

        subs = Clapeyron.split_model(model, [[1, 2, 4], [1, 3, 4]])
        @test subs[1].components == ["water", "sodium", "chloride"]
        @test subs[2].components == ["water", "potassium", "chloride"]
        @test water_activity(subs[1], 1e5, T, [1.0]) ≈ water_activity(model, salts, 1e5, T, [1.0, 0.0]) rtol = 1e-13

        model_fe = Pitzer(["water"], ["H+", "Fe3+", "Cl-"])
        miss_fe = pitzer_missing_mixing_parameters(model_fe)
        @test ("H+", "Fe3+") in miss_fe.theta
        @test ("H+", "Fe3+", "Cl-") in miss_fe.psi
    end

    @printline
end
end
