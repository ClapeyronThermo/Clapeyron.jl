using Test

include(joinpath(@__DIR__, "tp_flash_benchmark_cases.jl"))

cases = Dict(c.id => c for c in tpflash_benchmark_cases())

function _require_case(id::String)
    haskey(cases, id) || error("missing benchmark case: $id")
    return cases[id]
end

@testset "TP flash reference literals vs test_methods_api_flash literals" begin
    @testset "Scalar g⁺ literals" begin
        c = _require_case("vlle_pcsaft_water_cyclohexane_propane_298K_1bar")
        @test isapprox(c.g⁺, -6.759674475175065; rtol = 1e-6)

        c2 = _require_case("vle_pr_isobutane_nbutane_npenthane_nhexane_284p4K_1bar")
        @test isapprox(c2.g⁺, -6.618441125949686; rtol = 1e-6)

        c3 = _require_case("lle_pcsaft_water_cyclohexane_propane_298K_1bar")
        @test isapprox(c3.g⁺, -7.577270350886795; rtol = 1e-6)
    end

    @testset "Selected composition literals" begin
        c = _require_case("vle_454_case1_rr")
        @test isapprox(vec(c.x⁺[1, :]), [0.0023666624484214222, 0.9976333375515787, 0.0, 0.0]; rtol = 1e-6)

        c2 = _require_case("vle_454_case2_rr")
        @test isapprox(vec(c2.x⁺[1, :]), [5.306960867808201e-15, 0.010962897986743346, 0.6212133688148559, 0.3678237331983954]; rtol = 1e-6)
    end

    @testset "Selected β⁺ literals" begin
        # From `res2.fractions[1] ≈ 0.9991083897702745` in `test/test_methods_api_flash.jl`
        c = _require_case("vle_pr_isobutane_nbutane_npenthane_nhexane_282p3K_1bar")
        @test isapprox(c.β⁺[1], 0.9991083897702745; rtol = 1e-6)

        # From the #454 block: single-phase outcomes (one phase fraction is exactly zero)
        c3 = _require_case("vle_454_case3_rr")
        @test c3.β⁺ == [0.0, 1.0]

        c4 = _require_case("vle_454_case4_rr")
        @test c4.β⁺ == [0.0, 1.0]

        c5 = _require_case("vle_454_case5_rr")
        @test c5.β⁺ == [1.0, 0.0]

        c6 = _require_case("vle_454_case6_michelsen")
        @test c6.β⁺ == [0.0, 1.0]

        c7 = _require_case("vle_454_case7_michelsen")
        @test c7.β⁺ == [1.0, 0.0]
    end
end
