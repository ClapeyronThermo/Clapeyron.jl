using Test
using Clapeyron, ClapeyronHANNA

@testset "HANNA" begin
    model = HANNA(["isobutanol","water"];userlocations=(Mw = [74.1216, 18.01528], smiles = ["CC(C)CO", "O"]))

    T = 300.
    x1 = 0.5
    gEᵣ = Clapeyron.excess_gibbs_free_energy(model, NaN, T, [x1, 1-x1]) ./ Clapeyron.Rgas() ./ T
    γ = activity_coefficient(model, NaN, T, [x1, 1-x1])
    @test gEᵣ ≈ 0.5067839193591077 rtol = 1e-6
    @test log(γ[1]) ≈ 0.3337774981008488 rtol = 1e-6
    @test log(γ[2]) ≈ 0.6797903406173663 rtol = 1e-6
end