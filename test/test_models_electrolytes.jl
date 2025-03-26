@testset "Electrolyte models" begin
    @printline

    let T = 298.15, V = 1e-3, m = 1;

    @testset "ConstRSP" begin
        system = ConstRSP()
        Clapeyron.dielectric_constant(system, V, T, z) ≈ 78.38484961 rtol = 1e-6
    end

    @testset "ZuoFurst" begin
        system = ZuoFurst()
        Clapeyron.dielectric_constant(system, V, T, z) ≈ 78.30270731614397 rtol = 1e-6
    end

    @testset "LinMixRSP" begin
        system = LinMixRSP(["water"],["sodium","chloride"])
        Clapeyron.dielectric_constant(system, V, T, z) ≈ 75.9364209972588 rtol = 1e-6
    end

    @testset "Schreckenberg" begin
        system = Schreckenberg(["water"],["sodium","chloride"])
        @test Clapeyron.dielectric_constant(system, 1.8e-5 , T, z) ≈ 76.05272460268088 rtol = 1e-6
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

    @testset "SAFTgammaEMie" begin
        system = SAFTgammaEMie(["water"],["sodium","acetate"])
        salts = [("sodium acetate",["sodium"=>1,"acetate"=>1])]
        z = molality_to_composition(system,salts,m)
        #@test Clapeyron.a_res(system, V, T, z) ≈ -5.090195753415607 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -5.400523595315808 rtol = 1e-6
        
        model = SAFTgammaEMie(["water"],["calcium","chloride"]) #issue 349
        @test model isa EoSModel
    end
    @printline
end
end
