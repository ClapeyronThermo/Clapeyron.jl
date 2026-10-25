@testset verbose = true "Potential models" begin
    z1 = Clapeyron.SA[1.0]
    @testset "PeTSLJ" begin
        system = PeTSLJ()
        pure_generic_test(system)
    end

    @testset "Pohl2023Mie6" begin
        system = Pohl2023Mie6(12.0)
        pure_generic_test(system)
    end

    @testset "Thol2016LJ" begin
        system = Thol2016LJ()
        pure_generic_test(system)
    end

    @testset "KolafaNezbedaLJ" begin
        system = KolafaNezbedaLJ()
        pure_generic_test(system)
    end

    @testset "Johnson1993LJ" begin
        system = Johnson1993LJ()
        pure_generic_test(system)
    end
    
    #=
    @testset "Fe-ANN-Mie" begin
        system1 = FeANNMie(12.0,6.0)
        system2 = FeANNMie(12.0,6.0,solid = true)
        @test Clapeyron.a_potential(system1,2.0,1.0) ≈ -2.010945465814089 rtol = 1e-6
        @test Clapeyron.a_potential(system2,2.0,1.0) ≈ -2.0329590261282533 rtol = 1e-6
        test_scales(system,1.0)
    end=#

    @testset "LJRef" begin
        system = LJRef(["methane"],ljmodel = KolafaNezbedaLJ())
        pure_generic_test(system)
    end
end #testset
