@testset "Mass fraction inputs" begin
    model = PR(["methane", "ethane"])
    mass_fractions = MassFractions([0.25, 0.75],molar_basis = false)
    mass_amounts = [0.25u"kg", 0.75u"kg"]

    @test Clapeyron.uzstrip(model, mass_fractions) ≈ Clapeyron.uzstrip(model, mass_amounts)
    @test volume(model, 1u"bar", 300u"K", mass_fractions) ≈ volume(model, 1u"bar", 300u"K", mass_amounts)
end
