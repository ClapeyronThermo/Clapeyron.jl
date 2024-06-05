@testset "Electrolyte models" begin
    @printline
    
    let T = 298.15, V = 1e-3, m = 1;
    @testset "ePCSAFT" begin
        system = ePCSAFT(["water"],["sodium","chloride"])
        salts = [("sodium chloride",["sodium"=>1,"chloride"=>1])]
        z = molality_to_composition(system,salts,m)
        @test Clapeyron.a_res(system, V, T, z) ≈ -1.1094732499161208 rtol = 1e-6
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
        @test Clapeyron.a_res(system, V, T, z) ≈ -5.090195753415607 rtol = 1e-6
    end
    @printline
end
end
