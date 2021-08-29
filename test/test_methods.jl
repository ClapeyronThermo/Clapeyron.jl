
@testset "SAFT methods, single components" begin
    system = PCSAFT(["ethanol"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 5.907908736304141e-5 rtol = 1e-6 
        @test Clapeyron.volume(system, p, T;phase=:v) ≈ 0.020427920501436134 rtol = 1e-6 
        @test Clapeyron.volume(system, p, T;threaded=:false) ≈ 5.907908736304141e-5 rtol = 1e-6 
        @test Clapeyron.pip(system, 5.907908736304141e-5, T, [1.]) ≈ 6.857076349623449 rtol = 1e-6
        @test Clapeyron.compressibility_factor(system, p, T) ≈ 0.002383223535444557 rtol = 1e-6
        @test Clapeyron.pressure(system, 5.907908736304141e-5, T) ≈ p rtol = 1e-6
        @test Clapeyron.entropy(system, p, T) ≈ -58.87118569239617 rtol = 1E-6
        @test Clapeyron.chemical_potential(system, p, T)[1] ≈ -18323.877542682934 rtol = 1E-6
        @test Clapeyron.internal_energy(system, p, T) ≈ 2.9796670214913704e7 rtol = 1E-6
        @test Clapeyron.enthalpy(system, p, T) ≈ -35876.32155687084 rtol = 1E-6
        @test Clapeyron.gibbs_free_energy(system, p, T) ≈ -18323.87754268292 rtol = 1E-6
        @test Clapeyron.helmholtz_free_energy(system, p, T) ≈ -18329.785451419295 rtol = 1E-6
        @test Clapeyron.isochoric_heat_capacity(system, p, T) ≈ 48.37961296309505 rtol = 1E-6
        @test Clapeyron.isobaric_heat_capacity(system, p, T) ≈ 66.45719988319257 rtol = 1E-6
        @test Clapeyron.isothermal_compressibility(system, p, T) ≈ 1.1521981407243432e-9 rtol = 1E-6
        @test Clapeyron.isentropic_compressibility(system, p, T) ≈ 8.387789464951438e-10 rtol = 1E-6
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 1236.4846683094133 rtol = 1E-6 #requires that the model has Mr
        @test Clapeyron.isobaric_expansivity(system, p, T) ≈ -0.0010874255138433413 rtol = 1E-6
        @test Clapeyron.joule_thomson_coefficient(system, p, T) ≈ -6.007581864883784e-7 rtol = 1E-6
        @test Clapeyron.second_virial_coefficient(system, T) ≈ -0.004883874325089262 rtol = 1E-6
        @test Clapeyron.inversion_temperature(system, 1.1e8) ≈ 824.4137805298458 rtol = 1E-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 7972.550405922014 rtol = 1E-6
        @test Clapeyron.enthalpy_vap(system, T) ≈ 41712.78521121877 rtol = 1E-6
        @test Clapeyron.acentric_factor(system) ≈ 0.5730309964718605 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 533.1324329774004 rtol = 1E-6 
    end
end


@testset "LJSAFT methods, single components" begin
    system = LJSAFT(["ethanol"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 5.8990680856791996e-5 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 7933.046853495474 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 533.4350720160273 rtol = 1E-6 
    end
end

@testset "BACKSAFT methods, single components" begin
    system = BACKSAFT(["decane"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 0.00015924416586849443 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 581.785675150425 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 618.8455740197799 rtol = 1E-6 
    end
end

@testset "CPA methods, single components" begin
    system = CPA(["ethanol"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 5.913050998953597e-5 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 7923.883649594267 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 539.218257256262 rtol = 1E-6 
    end
end

@testset "SAFT-γ Mie methods, single components" begin
    system = SAFTγMie(["ethanol"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 5.753982584153832e-5 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 7714.849872968086 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 521.959273608428 rtol = 1E-6 
    end
end

@testset "SAFT-VR Mie methods, single components" begin
    system = SAFTVRMie(["methanol"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 4.064466003321247e-5 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 16957.625653548406 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 524.1487001618932 rtol = 1E-6 
    end
end

@testset "sCKSAFT methods, single component" begin
    system = sCKSAFT(["ethane"])
    tc_test,pc_test,vc_test = (321.00584034360014, 6.206975436514129e6, 0.0001515067748592245)
    tc,pc,vc = Clapeyron.crit_pure(system)
    @test tc ≈ tc_test rtol = 1E-3
    @test pc ≈ pc_test rtol = 1E-3
    @test vc ≈ vc_test rtol = 1E-3
end


@testset "SAFT methods, multi-components" begin
    system = PCSAFT(["methanol","cyclohexane"])
    p = 1e5
    T = 313.15
    z = [0.5,0.5]
    z_LLE = [0.27,0.73]
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T, z) ≈ 7.779694485714412e-5 rtol = 1e-6 
        @test Clapeyron.speed_of_sound(system, p, T, z) ≈ 1087.0303138908864 rtol = 1E-6
        @test Clapeyron.activity_coefficient(system, p, T, z)[1] ≈ 1.794138454452822 rtol = 1E-6
        @test Clapeyron.fugacity_coefficient(system, p, T, z)[1] ≈ 0.5582931304564298 rtol = 1E-6
    end
    @testset "Equilibrium properties" begin
        @test Clapeyron.bubble_pressure(system,T,z)[1] ≈ 54532.249600937736 rtol = 1E-6
        @test Clapeyron.LLE_pressure(system,T,z_LLE)[1] ≈ 737971.7522006684 rtol = 1E-6
        @test Clapeyron.three_phase(system, T)[1] ≈ 54504.079665621306 rtol = 1E-6
        @test Clapeyron.crit_mix(system,z)[1] ≈ 518.0004062881115 rtol = 1E-6
    end
end

@testset "Cubic methods, single components" begin
    system = RK(["ethane"])
    p = 1e7
    p2 = 1e5
    T = 250.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 6.819297582048736e-5 rtol = 1e-6 
        @test Clapeyron.volume(system, p2, T) ≈ 0.020539807199804024 rtol = 1e-6
        @test Clapeyron.volume(system, p2, T;phase=:vapour) ≈ 0.020539807199804024 rtol = 1e-6  
        @test Clapeyron.volume(system, p2, T;phase=:liquid) ≈ 7.563111462588624e-5 rtol = 1e-6 
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 800.288303407983 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 1.409820798879772e6 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 305.31999999999994 rtol = 1E-6 
    end
end

@testset "Cubic methods, multi-components" begin
    system = RK(["ethane","undecane"])
    p = 1e7
    T = 298.15
    z = [0.5,0.5]
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T, z) ≈ 0.00017378014541520907 rtol = 1e-6 
        @test Clapeyron.speed_of_sound(system, p, T, z) ≈ 892.4941848133369 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.bubble_pressure(system, T, z)[1] ≈ 1.5760730143760687e6 rtol = 1E-6
        @test Clapeyron.crit_mix(system, z)[1] ≈ 575.622237585033 rtol = 1E-6 
    end
end

@testset "Activity methods, multi-components" begin
    system = Wilson(["methanol","benzene"])
    p = 1e7
    T = 298.15
    z = [0.5,0.5]
    @testset "VLE properties" begin
        @test Clapeyron.bubble_pressure(system, T, z)[1] ≈ 23758.647133460465 rtol = 1E-6
    end
end

@testset "GERG2008 methods, single components" begin
    system = GERG2008(["water"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 1.8067969591040684e-5 rtol = 1e-6 
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 1484.0034692716843 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 3184.83242429761 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 647.0960000000457 rtol = 1E-6 
    end
end

@testset "GERG2008 methods, multi-components" begin
    @testset "Bulk properties" begin
        model = Clapeyron.GERG2008(["nitrogen","methane","ethane","propane","butane","isobutane","pentane"])
        lng_composition = [0.93,92.1,4.64,1.7,0.42,0.32,0.09]
        lng_composition_molar_fractions = lng_composition ./sum(lng_composition)
        @test Clapeyron.molar_density(model,(380.5+101.3)*1000.0,-153.0+273.15,lng_composition_molar_fractions)/1000 ≈ 24.98 rtol = 1E-2
        @test Clapeyron.mass_density(model,(380.5+101.3)*1000.0,-153.0+273.15,lng_composition_molar_fractions) ≈ 440.73 rtol = 1E-2
        @test Clapeyron.molar_density(model,(380.5+101.3)u"kPa",-153.0u"°C",lng_composition_molar_fractions;output=u"mol/L") ≈ 24.98*u"mol/L"  rtol=1E-2
        @test Clapeyron.mass_density(model,(380.5+101.3)u"kPa",-153.0u"°C",lng_composition_molar_fractions;output=u"kg/m^3")  ≈ 440.73*u"kg/m^3" rtol=1E-2
    end
    @testset "VLE properties" begin
        system = GERG2008(["carbon dioxide","water"])
        T = 298.15
        z = [0.8,0.2]
        @test Clapeyron.bubble_pressure(system, T,z)[1] ≈ 5.853909891112583e6 rtol = 1E-6
    end
end

@testset "IAPWS95 methods" begin
    system = IAPWS95()
    p = 1e5
    T = 298.15
    T_v = 380.15
    T_c = 750.
    p_c = 250e5
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 1.8068623941501927e-5 rtol = 1e-6 
        @test Clapeyron.volume(system, p, T_v;phase=:vapour) ≈ 0.03116877990373624 rtol = 1e-6 
        @test Clapeyron.volume(system, p_c, T_c) ≈ 0.00018553711945962424 rtol = 1e-6 
        @test Clapeyron.volume(system, p_c, T_c;phase=:sc) ≈ 0.00018553711945962424 rtol = 1e-6 
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 1496.699163371358 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 3169.9293390134403 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 647.096 rtol = 1E-5 
    end
end

@testset "PropaneRef methods" begin
    system = PropaneRef()
    p = 1e5
    T = 230.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 7.577761282115866e-5 rtol = 1e-6 
        @test Clapeyron.volume(system, p, T;phase=:vapour) ≈ 0.018421882342664616 rtol = 1e-6 
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 1166.6704395959607 rtol = 1e-6 
    end
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T)[1] ≈ 97424.11102152328 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 369.8900089509652 rtol = 1E-6 
    end
end

@testset "SPUNG methods" begin
    system = SPUNG(["ethane"])
    p = 1e5
    T = 313.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 0.035641472902311774 rtol = 1e-6 
        @test Clapeyron.volume(system, p, T;phase=:vapour) ≈ 0.035641472902311774 rtol = 1e-6 
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 357.8705332163255 rtol = 1e-6 
    end

    T_sat = 250.15
    @testset "VLE properties" begin
        @test Clapeyron.sat_pure(system, T_sat)[1] ≈ 3.5120264571020138e6 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 270.27247485012657 rtol = 1E-6 
    end
end