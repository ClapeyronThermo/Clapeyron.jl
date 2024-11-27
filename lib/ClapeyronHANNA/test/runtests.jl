using Test
using Clapeyron, ClapeyronHANNA

@testset "HANNA" begin
    T = 300.
    x1 = 0.5
    
    # Test 1
    model_1 = HANNA(["isobutanol","water"])

    gEᵣ_1 = Clapeyron.excess_gibbs_free_energy(model_1, NaN, T, [x1, 1-x1]) ./ Clapeyron.Rgas() ./ T
    γ_1 = activity_coefficient(model_1, NaN, T, [x1, 1-x1])
    @test gEᵣ_1 ≈ 0.5331872382421362 rtol = 1e-6
    @test log(γ_1[1]) ≈ 0.3235024331433443 rtol = 1e-6
    @test log(γ_1[2]) ≈ 0.7428720433409278 rtol = 1e-6

    # Test 2
    model_2 = HANNA(["diethylzinc","chloroform"],userlocations=(;Mw=[123.5,119.37],canonicalsmiles=["[CH2-]C.[CH2-]C.[Zn+2]","ClC(Cl)Cl"]))

    gEᵣ_2 = Clapeyron.excess_gibbs_free_energy(model_2, NaN, T, [x1, 1-x1]) ./ Clapeyron.Rgas() ./ T
    γ_2 = activity_coefficient(model_2, NaN, T, [x1, 1-x1])
    @test gEᵣ_2 ≈ -0.06662785416779005 rtol = 1e-5
    @test log(γ_2[1]) ≈ -0.06497140384559216 rtol = 1e-5
    @test log(γ_2[2]) ≈ -0.06828430448998822 rtol = 1e-5
end
