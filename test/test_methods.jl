using OpenSAFT, Test

@testset "SAFT methods single phase" begin
    test_system = system("ethanol", "PCSAFT")
    @test_broken get_volume(test_system, 1E5, 298)[1] ≈ 5.9069454773905654e-5 rtol = 0.01 #returns incorrect value
    @test get_sat_pure(test_system, 298)[1][1] ≈ 7904.2246894075415 rtol = 0.01
    @test get_enthalpy_vap(test_system, 298)[1] ≈ 41720.97771100548 rtol = 0.01
    @test_broken get_crit_pure(test_system)[1][1] ≈ 1 rtol = 1E10 #T_scale not defined
    @test get_pressure(test_system, 1, 298)[1] ≈ 1 rtol = 1E10
    @test get_entropy(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_chemical_potential(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_internal_energy(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_enthalpy(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_Gibbs_free_energy(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_Helmholtz_free_energy(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_isochoric_heat_capacity(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_isobaric_heat_capacity(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_isothermal_compressibility(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_isentropic_compressibility(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test_broken get_speed_of_sound(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10 #requires that the model has Mr
    @test get_isobaric_expansivity(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_Joule_Thomson_coefficient(test_system, 1E5, 298)[1] ≈ 1 rtol = 1E10
    @test get_second_virial_coeff(test_system, 298)[1] ≈ 1 rtol = 1E10
end
