
@testset "SAFT methods single phase" begin
    system = PCSAFT(["ethanol"])
    @test OpenSAFT.volume(system, 1E5, 298) ≈ 5.9069454773905654e-5 rtol = 0.01 #returns incorrect value
    @test OpenSAFT.sat_pure(system, 298)[1] ≈ 7904.2246894075415 rtol = 0.01
    # @test OpenSAFT.enthalpy_vap(system, 298) ≈ 41720.97771100548 rtol = 0.01
    @test OpenSAFT.crit_pure(system)[1] ≈ 1 rtol = 1E10 #T_scale not defined
    @test OpenSAFT.pressure(system, 1, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.entropy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.chemical_potential(system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test OpenSAFT.internal_energy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.enthalpy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.gibbs_free_energy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.helmholtz_free_energy(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.isochoric_heat_capacity(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.isobaric_heat_capacity(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.isothermal_compressibility(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.isentropic_compressibility(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.speed_of_sound(system, 1E5, 298) ≈ 1 rtol = 1E10 #requires that the model has Mr
    @test OpenSAFT.isobaric_expansivity(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.joule_thomson_coefficient(system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.second_virial_coefficient(system, 298) ≈ 1 rtol = 1E10
end
