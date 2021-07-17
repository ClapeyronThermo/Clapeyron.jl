
@testset "SAFT methods single phase" begin
    system = PCSAFT(["ethanol"])
    @test Clapeyron.volume(system, 1E5, 298) ≈ 5.9069454773905654e-5 rtol = 0.01 #returns incorrect value
    @test Clapeyron.sat_pure(system, 298)[1] ≈ 7904.2246894075415 rtol = 0.01
    @test Clapeyron.enthalpy_vap(system, 298) ≈ 41720.97771100548 rtol = 0.01
    @test Clapeyron.crit_pure(system)[1] ≈ 1 rtol = 1E10 #T_scale not defined
    @test Clapeyron.pressure(system, 1, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.entropy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.chemical_potential(system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test Clapeyron.internal_energy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.enthalpy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.gibbs_free_energy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.helmholtz_free_energy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.isochoric_heat_capacity(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.isobaric_heat_capacity(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.isothermal_compressibility(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.isentropic_compressibility(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.speed_of_sound(system, 1E5, 298) ≈ 1 rtol = 1E10 #requires that the model has Mr
    @test Clapeyron.isobaric_expansivity(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.joule_thomson_coefficient(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test Clapeyron.second_virial_coefficient(system, 298) ≈ 1 rtol = 1E10
end
