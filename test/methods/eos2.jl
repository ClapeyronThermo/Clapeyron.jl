
@testset verbose = true "Potential models  - single component methods" begin
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

    @testset "LJRef" begin
        system = LJRef(["methane"],ljmodel = KolafaNezbedaLJ())
        pure_generic_test(system)
    end
end #testset
