@testset verbose = true "Potential models" begin
    z1 = Clapeyron.SA[1.0]
    @testset "PeTSLJ" begin
        system = PeTSLJ()
        @test Clapeyron.a_potential(system,2.0,1.0) ≈ -1.5219558045689825 rtol = 1e-6
        test_scales(system,1.0)
    end

    @testset "Pohl2023Mie6" begin
        system = Pohl2023Mie6(12.0)
        @test Clapeyron.a_potential(system,2.0,1.0) ≈ -2.0551233779172424 rtol = 1e-6
        test_scales(system,1.0)

        #suplementary information test values
        @test Clapeyron.pressure(Pohl2023Mie6(11.0),100,0.9) ≈ 0.008374117 rtol = 1e-6
        @test Clapeyron.pressure(Pohl2023Mie6(11.0),1/0.4,1.5) ≈ 0.226112187 rtol = 1e-6

        @test Clapeyron.pressure(Pohl2023Mie6(12.0),100,0.9) ≈ 0.008425582 rtol = 1e-6
        @test Clapeyron.pressure(Pohl2023Mie6(12.0),1/0.4,1.5) ≈ 0.283040888 rtol = 1e-6

        @test Clapeyron.pressure(Pohl2023Mie6(13.0),100,0.9) ≈ 0.008468062 rtol = 1e-6
        @test Clapeyron.pressure(Pohl2023Mie6(13.0),1/0.4,1.5) ≈ 0.330515553 rtol = 1e-6
    end

    @testset "Thol2016LJ" begin
        system = Thol2016LJ()
        @test Clapeyron.a_potential(system,2.0,1.0) ≈ -2.0501021202423115 rtol = 1e-6
        test_scales(system,1.0)
    end

    @testset "KolafaNezbedaLJ" begin
        system = KolafaNezbedaLJ()
        @test Clapeyron.a_potential(system,2.0,1.0) ≈ -2.052843069659619 rtol = 1e-6
        test_scales(system,1.0)
    end

    @testset "Johnson1993LJ" begin
        system = Johnson1993LJ()
        @test Clapeyron.a_potential(system,2.0,1.0) ≈ -2.039826585613582 rtol = 1e-6
        test_scales(system,1.0)
    end
    
    
    @testset "Fe-ANN-Mie" begin
        system1 = FeANNMie(12.0,6.0)
        system2 = FeANNMie(12.0,6.0,solid = true)
        @test Clapeyron.a_potential(system1,2.0,1.0) ≈ -2.010945465814089 rtol = 1e-6
        @test Clapeyron.a_potential(system2,2.0,1.0) ≈ -2.0329590261282533 rtol = 1e-6
        test_scales(system1,1.0)
    end

    @testset "LJRef" begin
        system = LJRef(["methane"])
        f,h = inverse_shape_factors(system,1.0)
        Tx = 1.051*f
        f2,h2 = shape_factors(system,f)
        @test f2 ≈ f
        Vx = h/0.673
        z1 = Clapeyron.SA[1.0]
        @test Clapeyron.a_res(system, Vx, Tx, z1) ≈ -2.244730279521925 rtol = 1e-6
    end
end #testset
