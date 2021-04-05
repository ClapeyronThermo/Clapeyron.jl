using OpenSAFT, Test

@testset "SAFT methods single phase" begin
    test_system = PCSAFT(["ethanol"])
    @test OpenSAFT.volume(test_system, 1E5, 298) ≈ 5.9069454773905654e-5 rtol = 0.01 #returns incorrect value
    @test OpenSAFT.sat_pure(test_system, 298)[1] ≈ 7904.2246894075415 rtol = 0.01
    @test OpenSAFT.enthalpy_vap(test_system, 298) ≈ 41720.97771100548 rtol = 0.01
    @test_broken OpenSAFT.crit_pure(test_system)[1] ≈ 1 rtol = 1E10 #T_scale not defined
    @test OpenSAFT.pressure(test_system, 1, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.entropy(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.chemical_potential(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test OpenSAFT.internal_energy(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.enthalpy(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.gibbs_free_energy(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.helmholtz_free_energy(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.isochoric_heat_capacity(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.isobaric_heat_capacity(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.isothermal_compressibility(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.isentropic_compressibility(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.speed_of_sound(test_system, 1E5, 298) ≈ 1 rtol = 1E10 #requires that the model has Mr
    @test OpenSAFT.isobaric_expansivity(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.joule_thomson_coefficient(test_system, 1E5, 298) ≈ 1 rtol = 1E10
    @test OpenSAFT.second_virial_coeff(test_system, 298) ≈ 1 rtol = 1E10
end
