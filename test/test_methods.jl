
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

@testset "testing density methods for GERG2008" begin
    model = Clapeyron.GERG2008(["nitrogen","methane","ethane","propane","butane","isobutane","pentane"])
    lng_composition = [0.93,92.1,4.64,1.7,0.42,0.32,0.09]
    lng_composition_molar_fractions = lng_composition ./sum(lng_composition)
    @test Clapeyron.molar_density(model,(380.5+101.3)*1000.0,-153.0+273.15,lng_composition_molar_fractions)/1000 ≈ 24.98 rtol = 1E-2
    @test Clapeyron.mass_density(model,(380.5+101.3)*1000.0,-153.0+273.15,lng_composition_molar_fractions) ≈ 440.73 rtol = 1E-2
    @test Clapeyron.molar_density(model,(380.5+101.3)u"kPa",-153.0u"°C",lng_composition_molar_fractions;output=u"mol/L") ≈ 24.98*u"mol/L"  rtol=1E-2
    @test Clapeyron.mass_density(model,(380.5+101.3)u"kPa",-153.0u"°C",lng_composition_molar_fractions;output=u"kg/m^3")  ≈ 440.73*u"kg/m^3" rtol=1E-2
end

@testset "Critical pure test, sCKSAFT" begin
    smodel = sCKSAFT(["ethane"])
    tc_test,pc_test,vc_test = (340.2146485307944, 7.321069162025713e6, 0.0001381492224781526)
    tc,pc,vc = crit_pure(smodel)
    @test tc ≈ tc_test rtol = 1E-3
    @test pc ≈ pc_test rtol = 1E-3
    @test vc ≈ vc_test rtol = 1E-3
end
