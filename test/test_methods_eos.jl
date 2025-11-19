using Clapeyron, Test, Unitful

@testset "pharmaPCSAFT, single components" begin
    system = pharmaPCSAFT(["water"])
    v1 = Clapeyron.saturation_pressure(system, 280.15)[2]
    v2 = Clapeyron.saturation_pressure(system, 278.15)[2]
    v3 = Clapeyron.saturation_pressure(system, 275.15)[2]

    @test v1 ≈ 1.8022929328333385e-5  rtol = 1E-6
    @test v2 ≈ 1.8022662044726256e-5 rtol = 1E-6
    @test v3 ≈ 1.802442451376152e-5 rtol = 1E-6
    #density maxima of water
    @test v2 < v1
    @test v2 < v3

    #issue 377
    mod_phsft = pharmaPCSAFT(["water", "hydrogen"])
    p_c = 3e6
    T = 343.15
    y = [1.955278169111263e-5, 0.9999804472183089]
    @test Clapeyron.volume(mod_phsft,p_c,T,y,phase = :v) ≈ 0.0009700986016167609 rtol = 1E-6
end

@testset "LJSAFT methods, single components" begin
    system = LJSAFT(["ethanol"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 5.8990680856791996e-5 rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 7933.046853495474 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 533.4350720160273 rtol = 1E-6
    end
end

@testset "softSAFT methods, single components" begin
    system = softSAFT(["ethanol"])
    solid_system = solidsoftSAFT("octane")
    p = 1e5
    T = 273.15 + 78.24
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T,phase=:v) ≈ 0.027368884099868623 rtol = 1e-6
        #volume(SAFTgammaMie(["ethanol"]),p,T,phase =:l)  =6.120507339375205e-5
        @test Clapeyron.volume(system, p, T,phase=:l) ≈ 6.245903786961202e-5 rtol = 1e-6
        @test Clapeyron.volume(solid_system,2.3e9,298.15,phase = :s) ≈ 9.961905037894007e-5 rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 101341.9709136089 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 540.1347889779657 rtol = 1E-6
    end

end

@testset "BACKSAFT methods, single components" begin
    system = BACKSAFT(["decane"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        #0.0001950647173402879 with SAFTgammaMie
        @test Clapeyron.volume(system, p, T) ≈ 0.00019299766073894634 rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 167.8313793818096 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 618.8455740197799 rtol = 1E-6
    end
end

@testset "CPA methods, single components" begin
    system = CPA(["ethanol"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 5.913050998953597e-5 rtol = 1e-6
        @test volume(CPA("water"), 1e5u"Pa", 303.15u"K") ≈ 1.7915123921401366e-5u"m^3" rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 7923.883649594267 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 539.218257256262 rtol = 1E-6
    end
end

@testset "SAFT-γ Mie methods, single components" begin
    system = SAFTγMie(["ethanol"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 5.753982584153832e-5 rtol = 1e-6
        @test Clapeyron.molecular_weight(system)*1000 ≈ 46.065
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 7714.8637084302 rtol = 1E-5
        @test Clapeyron.crit_pure(system)[1] ≈ 521.963002384691 rtol = 1E-5
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
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 16957.59261579083 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 519.5443602179253  rtol = 1E-5
    end
end

@testset "SAFT-VRQ Mie methods, single component" begin
    system = SAFTVRQMie(["helium"])
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, 4)[1] ≈ 56761.2986265459 rtol = 1E-6
    end
end

@testset "sCKSAFT methods, single component" begin
    system = sCKSAFT(["ethane"])
    tc_test,pc_test,vc_test = (321.00584034360014, 6.206975436514129e6, 0.0001515067748592245)
    tc,pc,vc = Clapeyron.crit_pure(system)
    @test tc ≈ tc_test rtol = 1E-3
    @test pc ≈ pc_test rtol = 1E-3
    @test vc ≈ vc_test rtol = 1E-3
    system2 = CKSAFT("methanol")
    @test crit_pure(system2)[1] ≈ 544.4777700204786
end

@testset "ideal model parsing to JSON" begin
    comp = ["hexane"]
    idealmodels = []
    push!(idealmodels,BasicIdeal(comp))
    push!(idealmodels,MonomerIdeal(comp))
    push!(idealmodels,JobackIdeal(comp))
    push!(idealmodels,WalkerIdeal(comp))
    push!(idealmodels,AlyLeeIdeal(comp))

    T0 = 400.15
    V0 = 0.03
    for mi in idealmodels
        mxd = XiangDeiters(comp, idealmodel = mi)
        id_mxd = Clapeyron.idealmodel(mxd)
       
        #parsed and reconstituted idealmodel
        cp1 = Clapeyron.VT_isobaric_heat_capacity(id_mxd,V0,T0)
        #original
        a1 = Clapeyron.a_ideal(id_mxd,V0,T0,Clapeyron.SA[1.0])
        a2 = Clapeyron.a_ideal(mi,V0,T0,Clapeyron.SA[1.0])
        cp2 = Clapeyron.VT_isobaric_heat_capacity(mi,V0,T0)
        if mi isa AlyLeeIdeal
            @test_broken a1 ≈ a2 rtol = 1e-6
            @test_broken cp1 ≈ cp2 rtol = 1e-6
        else
            @test a1 ≈ a2 rtol = 1e-6
            @test cp1 ≈ cp2 rtol = 1e-6
        end
    end
end

@testset "SAFT-VRQ Mie methods, multicomponent" begin
    system = SAFTVRQMie(["hydrogen","neon"])
    T = -125 + 273.15
    #brewer 1969 data for H2-Ne
    #Texp = [50,25,0,-25,-50,-75,-100,-125] .+ 273.15
    #B12exp = [14.81,14.23,13.69,13.03,12.29,11.23,9.98,8.20]
    #there is a difference between brewer 1969 data and the exact value, but for some reason, their plots use a very thick linewidth...
    @test Clapeyron.equivol_cross_second_virial(system, T)*1e6 ≈ 8.09 rtol = 1E-1
end

@testset "RK, single component" begin
    system = RK(["ethane"])
    p = 1e7
    p2 = 1e5
    T = 250.15
    @testset "Bulk properties" begin
        #Check that we actually dispatch to volume_impl, if they vary by anything, we are using the default volume solver
        @test Clapeyron.volume(system, p, T) == Clapeyron.volume_impl(system, p, T)
        @test Clapeyron.volume_impl(system, p, T) ≈ 6.819297582048736e-5 rtol = 1e-6
        @test Clapeyron.volume_impl(system, p2, T) ≈ 0.020539807199804024 rtol = 1e-6
        @test Clapeyron.volume_impl(system, p2, T,[1.0], :vapour) ≈ 0.020539807199804024 rtol = 1e-6
        @test Clapeyron.volume_impl(system, p2, T, [1.0], :liquid) ≈ 7.563111462588624e-5 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 800.288303407983 rtol = 1e-6
    end
    @testset "VLE properties" begin
        psat,_,_ = Clapeyron.saturation_pressure(system, T)
        @test psat ≈ 1.409820798879772e6 rtol = 1E-6
        @test Clapeyron.saturation_pressure(system, T,SuperAncSaturation())[1]  ≈ psat rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 305.31999999999994 rtol = 1E-6
        @test Clapeyron.wilson_k_values(system,p,T) ≈ [0.13840091523637849]  rtol = 1E-6
    end
end
GC.gc()
@testset "Patel-Teja, single component" begin
    system = PatelTeja(["water"])
    p = 1e5
    T = 550.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 0.04557254632681239 rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 4.51634223156497e6 rtol = 1E-6
        Tc,Pc,Vc = Clapeyron.crit_pure(system)
        @test Tc == system.params.Tc.values[1]
        @test Pc == system.params.Pc.values[1]
        @test pressure(system,Vc,Tc) ≈ Pc
    end
end

@testset "Patel-Teja-Valderrama, single component" begin
    system = PTV(["water"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 1.9221342043684064e-5 rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 2397.1315826665273 rtol = 1E-6
        Tc,Pc,Vc = Clapeyron.crit_pure(system)
        @test Tc == system.params.Tc.values[1]
        @test pressure(system,Vc,Tc) ≈ Pc
    end
end

@testset "KU, single component" begin
    system = KU(["water"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 2.0614093101238483e-5 rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 2574.7348636211996 rtol = 1E-6
    end
end

@testset "RKPR, single component" begin
    system = RKPR(["methane"])
    vc_vol = volume(system,system.params.Pc[1],system.params.Tc[1]) #vc calculated via cubic_poly
    Tc,Pc,Vc = crit_pure(system) #vc calculated via cubic_pure_zc
    @test vc_vol ≈ Vc rtol = 1e-4
    @test vc_vol/system.params.Vc[1] ≈ 1.168 rtol = 1e-4 #if Zc_exp < 0.29, this should hold, by definition
    @test Vc/system.params.Vc[1] ≈ 1.168 rtol = 1e-4 #if Zc_exp < 0.29, this should hold, by definition
end

@testset "EPPR78, single component" begin
    system = EPPR78(["carbon dioxide"])
    T = 400u"K"
    @test Clapeyron.volume(system, 3311.0u"bar", T) ≈ 3.363141761376883e-5u"m^3"
    @test Clapeyron.molar_density(system, 3363.1u"bar", T) ≈ 29810.09484964839u"mol*m^-3"
end

@testset "Cubic methods, multi-components" begin
    system = RK(["ethane","undecane"])
    system2 = tcPR(["benzene","toluene","nitrogen"])
    system3 = cPR(["butane","toluene"],idealmodel = ReidIdeal)
    p = 1e7
    T = 298.15
    z = [0.5,0.5]
    p2 = 1.5*101325
    T2 = 350
    z2 = [0.001,0.001,0.001]

    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T, z) ≈ 0.00017378014541520907 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T, z) ≈ 892.4941848133369 rtol = 1e-6
        @test Clapeyron.volume(system2, p2, T2, z2, phase = :l) ≈ 2.851643999862116e-7 rtol = 1e-6
        @test Clapeyron.volume(system2, p2, T2, z2 ./ sum(z2), phase = :l) ≈ 2.851643999862116e-7/sum(z2) rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.bubble_pressure(system, T, z)[1] ≈ 1.5760730143760687e6 rtol = 1E-6
        @test Clapeyron.crit_mix(system, z)[1] ≈ 575.622237585033 rtol = 1E-6
        @test Clapeyron.mechanical_critical_point(system,z)[1] ≈ 483.08783138464617 rtol = 1E-6
        @test Clapeyron.spinodal_maximum(system,z)[1] ≈ 578.7715447554762 rtol = 1E-6
        srksystem  = SRK(["ethane","undecane"])
        @test Clapeyron.wilson_k_values(srksystem,p,T) ≈ [0.4208525854463047, 1.6171551943938252e-5] rtol = 1E-6
        #test scaling of crit_mix
        cm1 = crit_mix(system3,[0.5,0.5])
        cm2 = crit_mix(system3,[1.0,1.0])
        @test cm1[1] ≈ cm2[1]
        @test 2*cm1[3] ≈ cm2[3]
    end
end
GC.gc()
@testset "Activity methods, pure components" begin
    if hasfield(Wilson,:puremodel)
        system = Wilson(["methanol"])
    else
        system = CompositeModel(["methanol"],liquid = Wilson,fluid = PR)
    end
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 4.7367867309516085e-5 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 2136.222735675237 rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.crit_pure(system)[1] ≈ 512.6400000000001 rtol = 1E-6
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 15525.93630485447 rtol = 1E-6
    end
end

@testset "Activity methods, multi-components" begin
    com = CompositeModel(["water","methanol"],liquid = DIPPR105Liquid,saturation = DIPPR101Sat,gas = PR)
    
    system = Wilson(["methanol","benzene"])
    comp_system = CompositeModel(["methanol","benzene"]; fluid = PR, liquid = Wilson,reference_state = :ashrae)


    if hasfield(Wilson,:puremodel)
        system2 = Wilson(["water","methanol"],puremodel = com)
    else
        system2 = CompositeModel(["water","methanol"],liquid = Wilson, fluid = com)
    end

    if hasfield(UNIFAC,:puremodel)
        system3 = UNIFAC(["octane","heptane"],puremodel = LeeKeslerSat)
    else
        system3 = CompositeModel(["octane","heptane"],liquid = UNIFAC,fluid = LeeKeslerSat)
    end

    com1 = split_model(com)[1]
    p = 1e5
    T = 298.15
    T2 = 320.15
    z = [0.5,0.5]
    z_bulk = [0.2,0.8]
    T3 = 300.15
    z3 = [0.9,0.1]
    @testset "Bulk properties" begin
        @test crit_pure(com1)[1] ≈ 647.13
        @test Clapeyron.volume(system, p, T, z_bulk) ≈ 7.967897222918716e-5 rtol = 1e-6
        @test Clapeyron.volume(comp_system, p, T, z_bulk) ≈ 7.967897222918716e-5 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T, z_bulk) ≈ 1551.9683977722198 rtol = 1e-6
        @test Clapeyron.speed_of_sound(comp_system, p, T, z_bulk) ≈ 1551.9683977722198 rtol = 1e-6
        @test Clapeyron.mixing(system, p, T, z_bulk, Clapeyron.gibbs_free_energy) ≈ -356.86007792929263 rtol = 1e-6
        @test Clapeyron.mixing(system, p, T, z_bulk, Clapeyron.enthalpy) ≈ 519.0920708672975 rtol = 1e-6
        #test that we are actually considering the reference state, even in the vapour phase.
        @test enthalpy(comp_system,p,T,z_bulk,phase = :v) - enthalpy(system,p,T,z_bulk,phase = :v) ≈ sum(reference_state(comp_system).a0 .* z_bulk) rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.gibbs_solvation(system, T) ≈ -24707.145697543132 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system, T, z)[1] ≈ 23758.58099358788 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system, T, z,ActivityBubblePressure(gas_fug = true,poynting = true))[1] ≈ 23839.554959977086
        #@test Clapeyron.bubble_pressure(system, T, z,ActivityBubblePressure(gas_fug = true,poynting = false))[1] ≈ 23833.324475723246
        @test Clapeyron.bubble_pressure(system3,T3,z3)[1] ≈ 2460.897944633704 rtol = 1E-6
        @test Clapeyron.bubble_temperature(system,23758.58099358788, z)[1] ≈ T rtol = 1E-6

        @test Clapeyron.dew_pressure(system2, T2, z)[1] ≈ 19386.939256733036 rtol = 1E-6
        @test Clapeyron.dew_pressure(system2, T2, z,ActivityDewPressure(gas_fug = true,poynting = true))[1] ≈ 19393.924550078184 rtol = 1e-6
        #@test Clapeyron.dew_pressure(system2, T2, z,ActivityDewPressure(gas_fug = true,poynting = false))[1] ≈ 19393.76058757084 rtol = 1e-6
        #@test Clapeyron.dew_temperature(system2, 19386.939256733036, z)[1]  ≈ T2 rtol = 1E-6
    end

    @testset "LLE" begin
        model3 = NRTL(["methanol","hexane"])
        x1,x2 = Clapeyron.LLE(model3,290.0)
        @test x1[1] ≈ 0.15878439462531743 rtol = 1E-6
    end
end
GC.gc()
@testset "GERG2008 methods, single components" begin
    system = GERG2008(["water"])
    met = GERG2008(["methane"])
    p = 1e5
    T = 298.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 1.8067969591040684e-5 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 1484.0034692716843 rtol = 1e-6
        #EOS-LNG, table 15
        V1,T1 = 1/27406.6102,100.0
        @test Clapeyron.pressure(met,V1,T1)  ≈ 1.0e6 rtol = 2e-6
        @test Clapeyron.VT_speed_of_sound(met,V1,T1) ≈ 1464.5158 rtol = 1e-6
        @test Clapeyron.pressure(met,1/28000,140) ≈ 86.944725e6  rtol = 2e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 3184.83242429761 rtol = 1E-6
        @test Clapeyron.saturation_pressure(system, T, IsoFugacitySaturation())[1] ≈ 3184.83242429761 rtol = 1E-6
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
    
        #test found in #371
        model2 = GERG2008(["carbon dioxide","nitrogen","water"])
        @test mass_density(model2,64.0e5,30+273.15,[0.4975080785711593, 0.0049838428576813995, 0.4975080785711593],phase = :l) ≈ 835.3971524715569 rtol = 1e-6
    
        #test found in #395:

        p395 = range(log(3.12e6),log(1e8),1000)
        model395 = GERG2008(["carbon dioxide","nitrogen"])
        z395 = Ref([0.95,0.05])
        v395 = volume.(model395,exp.(p395),250.0,z395,phase = :v)
        @test count(isnan,v395) == 999

    end
    @testset "VLE properties" begin
        system = GERG2008(["carbon dioxide","water"])
        T = 298.15
        z = [0.8,0.2]
        @test Clapeyron.bubble_pressure(system, T,z)[1] ≈ 5.853909891112583e6 rtol = 1E-5
    end
end

@testset "EOS-LNG methods, multi-components" begin
    @testset "Bulk properties" begin
        #EOS-LNG paper, table 16
        system = EOS_LNG(["methane","isobutane"])
        z = [0.6,0.4]
        T1,V1 = 160.0,1/17241.868
        T2,V2 = 350.0,1/100

        @test Clapeyron.VT_speed_of_sound(system,V1,T1,z) ≈ 1331.9880 rtol = 1e-6
        @test Clapeyron.VT_speed_of_sound(system,V2,T2,z) ≈ 314.72845 rtol = 1e-6
        @test Clapeyron.volume(system,5e6,T1,z) ≈ V1
        @test Clapeyron.pressure(system,V2,T2,z) ≈ 0.28707693e6 rtol = 1e6
    end
end

@testset "IAPWS95 methods" begin
    system = IAPWS95()
    p = 1e5
    T = 298.15
    T_v = 380.15
    T_c = 750.
    p_c = 250e5
    mw = Clapeyron.molecular_weight(system)
    @testset "Bulk properties" begin
        #IAPWS-2018, table 7
        @test mass_density(system,0.992418352e5,300.0) ≈ 996.556 rtol = 1e-6
        @test mass_density(system,0.200022515e8,300.0) ≈ 1005.308 rtol = 1e-6
        @test mass_density(system,0.700004704e9,300.0) ≈ 1188.202 rtol = 1e-6
        test_volume(system,0.992418352e5,300.0)
        test_volume(system,0.200022515e8,300.0)
        test_volume(system,0.700004704e9,300.0)
        @test entropy(system,0.992418352e5,300.0) ≈ mw*393.062643 rtol = 1e-6
        @test entropy(system,0.200022515e8,300.0) ≈ mw*387.405401 rtol = 1e-6
        @test entropy(system,0.700004704e9,300.0) ≈ mw*132.609616 rtol = 1e-6
        @test speed_of_sound(system,0.992418352e5,300.0) ≈ 1501.51914 rtol = 1e-6
        @test speed_of_sound(system,0.200022515e8,300.0) ≈ 1534.92501 rtol = 1e-6
        @test speed_of_sound(system,0.700004704e9,300.0) ≈ 2443.57992 rtol = 1e-6
        #below triple point/ liquid
        test_volume(system,1e6,265.0,phase = :l)
        test_volume(system,1e8,265.0,phase = :l)

        @test Clapeyron.molecular_weight(Clapeyron.idealmodel(system)) == mw
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 3169.9293388718283  rtol = 1E-6
        @test Clapeyron.saturation_pressure(system, T, IsoFugacitySaturation())[1] ≈ 3169.9293388718283 rtol = 1E-6
        #saturation temperature tests are noisy
        @test Clapeyron.saturation_temperature(system,3169.9293390134403)[1] ≈ 298.1499999999789 rtol = 1E-6
        tc,pc,vc = Clapeyron.crit_pure(system)
        @test tc ≈ 647.096 rtol = 1E-5
        v2 = volume(system,pc,tc)
        @test pressure(system,v2,tc) ≈ pc rtol = 1E-6
    end
end

@testset "PropaneRef methods" begin
    system = PropaneRef()
    p = 1e5
    T = 230.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 7.577761282115866e-5 rtol = 1e-6
        #PropsSI("D","P",1e4,"T",230.15,"propane") == 0.23126803007122876
        @test Clapeyron.mass_density(system, 1e4, T) ≈ 0.23126803007122876 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 1166.6704395959607 rtol = 1e-6
    end
    @testset "VLE properties" begin
        ps = 97424.11102296013 #PropsSI("P","T",T,"Q",1,"propane")
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ ps rtol = 1E-6
        @test Clapeyron.saturation_pressure(system, T, IsoFugacitySaturation())[1] ≈ ps rtol = 1E-6
        @test Clapeyron.saturation_temperature(system,ps)[1] ≈ T  rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 369.8900089509652 rtol = 1E-6
    end
end

@testset "Ammonia2023 methods" begin
    system = Ammonia2023()
    #test for pressures here. we exclude the 0 density (obtained from paper's SI)
    ps = [
       27.63275
       11.25287
       0.4270652
       5.96888
       58.52754]
    Tρ = [(320.0,35.0),(405.0,16.5),(275.0,0.2),(520.0,1.5),(620.0,14.0)]
    for i in eachindex(ps)
        T,rho = Tρ[i]
        @test Clapeyron.pressure(system,0.001/rho,T)*1e-6 ≈ ps[i] rtol = 1e-6
    end

end
GC.gc()
@testset "Helmholtz + Activity" begin
    model = HelmAct(["water","ethanol"])
    p = 12666.0
    x1 = Clapeyron.FractionVector(0.00350)
    test_scales(model) #
    @test bubble_temperature(model,p,x1)[4][1] ≈ 0.00198 rtol = 1e-2
end

@testset "SingleFluid - CoolProp" begin
    #methanol, uses double exponential term
    #before, it used the association term, but now no model uses it
    @test saturation_pressure(SingleFluid("methanol"),300.15)[1] ≈ PropsSI("P","T",300.15,"Q",1.,"methanol") rtol = 1e-6
    
    r134 = SingleFluid("r134a")
    r1342 = MultiFluid("r134a")
    @test Clapeyron.eos(r134,0.03,373.15,Clapeyron.SA[1.0]) ≈ PropsSI("HELMHOLTZMOLAR","Dmolar",1/0.03,"T",373.15,"R134a")
    @test Clapeyron.eos(r1342,0.03,373.15,Clapeyron.SA[1.0]) ≈ PropsSI("HELMHOLTZMOLAR","Dmolar",1/0.03,"T",373.15,"R134a")
    @test Clapeyron.a_res(r134,0.03,373.15,Clapeyron.SA[1.0]) ≈ PropsSI("ALPHAR","Dmolar",1/0.03,"T",373.15,"R134a")
    @test Clapeyron.a_res(r1342,0.03,373.15,Clapeyron.SA[1.0]) ≈ PropsSI("ALPHAR","Dmolar",1/0.03,"T",373.15,"R134a")

    #tests send via email

    fluid1 = SingleFluid("n-Undecane")
    test_volume(fluid1,1e-2*fluid1.properties.Pc,0.38*fluid1.properties.Tc)
    test_volume(fluid1,3e2*fluid1.properties.Pc,0.38*fluid1.properties.Tc)
    test_volume(fluid1,3e2*fluid1.properties.Pc,1.1*fluid1.properties.Tc)

    fluid2 = SingleFluid("n-Butane")
    test_volume(fluid2,1e-2*fluid2.properties.Pc,0.3*fluid2.properties.Tc)
    test_volume(fluid2,30*fluid2.properties.Pc,0.3*fluid2.properties.Tc)

    fluid3 = SingleFluid("water")
    test_volume(fluid3,1e-2*fluid3.properties.Pc,0.4*fluid3.properties.Tc)
    test_volume(fluid3,40*fluid3.properties.Pc,3.2*fluid3.properties.Tc)

    fluid4 = SingleFluid("MethylOleate")
    test_volume(fluid4,1e-2*fluid4.properties.Pc,0.3*fluid4.properties.Tc)
    test_volume(fluid4,4e1*fluid4.properties.Pc,0.3*fluid4.properties.Tc)
    test_volume(fluid4,4e1*fluid4.properties.Pc,1.3*fluid4.properties.Tc)

    fluid5 = SingleFluid("MD3M")
    test_volume(fluid5,1e-2*fluid5.properties.Pc,0.3*fluid5.properties.Tc)
    test_volume(fluid5,2e2*fluid5.properties.Pc,0.3*fluid5.properties.Tc)
    test_volume(fluid5,2e2*fluid5.properties.Pc,1.1*fluid5.properties.Tc)

    fluid6 = SingleFluid("Toluene")
    test_volume(fluid6,1e-2*fluid6.properties.Pc,0.25*fluid6.properties.Tc)
    test_volume(fluid6,2e2*fluid6.properties.Pc,0.25*fluid6.properties.Tc)
    test_volume(fluid6,2e2*fluid6.properties.Pc,1.2*fluid6.properties.Tc)

    #CoolProp fluid predicting negative fundamental derivative of gas dynamics
    #10.1021/acs.iecr.9b00608, figure 17
    model = SingleFluid("MD4M")
    TΓmin = 647.72
    _,_,vv = saturation_pressure(model,TΓmin)
    Γmin = Clapeyron.VT_fundamental_derivative_of_gas_dynamics.(model,vv,TΓmin)
    @test Γmin ≈ -0.2825376983518102 rtol = 1e-6

    #376
    p1 = 1e5 .* (0.1:0.1:420)
    T_376 = (310.95,477.95,644.15)
    px = first.(saturation_pressure.(fluid3,T_376))
    p1 = px[1]:1e4:420e5
    p2 = px[2]:1e4:420e5
    p3 = px[3]:1e4:420e5
    v_T37 = volume.(fluid3,p1,T_376[1],phase = :l)
    v_T202 = volume.(fluid3,p2,T_376[2],phase = :l)
    v_T371 = volume.(fluid3,p3,T_376[3],phase = :l)
    @test iszero(count(isnan,v_T37))
    @test iszero(count(isnan,v_T202))
    @test iszero(count(isnan,v_T371))

    #pseudo pure
    pseudo_pure = EmpiricPseudoPure("R410A")
    @test bubble_pressure(pseudo_pure,220.0,[1.0])[1] ≈ PropsSI("P","T",220.0,"Q",0.,"R410A") rtol = 1e-6
    @test dew_pressure(pseudo_pure,220.0,[1.0])[1] ≈ PropsSI("P","T",220.0,"Q",1.,"R410A") rtol = 1e-6
    @test bubble_temperature(pseudo_pure,1e5,[1.0])[1] ≈ PropsSI("T","P",1e5,"Q",0.,"R410A") rtol = 1e-6
    @test dew_temperature(pseudo_pure,1e5,[1.0])[1] ≈ PropsSI("T","P",1e5,"Q",1.,"R410A") rtol = 1e-6

end

@testset "LKP methods" begin
    system = LKP("propane", idealmodel = AlyLeeIdeal)
    system_mod = LKPmod("squalane",userlocations = (Tc = 810,Pc = 0.728e6,acentricfactor = 1.075,Mw = 1.0))

    p = 1e5
    T = 230.15
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T, phase = :l) ≈ 7.865195401331961e-5 rtol = 1e-6
        @test Clapeyron.volume(system, p, T, phase = :v) ≈ 0.018388861273788176 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T, phase = :l) ≈ 1167.2461897307874 rtol = 1e-6
        @test Clapeyron.molar_density(system_mod,0.0,298.15,phase =:l) ≈ 1721.2987626107251 rtol = 1e-6 #0.1007/s10765-024-03360-0, Figure 4
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 105419.26772976149 rtol = 1E-6
        #saturation temperature tests are noisy
        @test Clapeyron.saturation_temperature(system,105419.26772976149)[1] ≈ T  rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 369.83 rtol = 1E-6
    end
end

@testset "LJRef methods" begin
    system = LJRef(["methane"])
    T = 1.051*Clapeyron.T_scale(system)
    p = 0.035*Clapeyron.p_scale(system)
    v = Clapeyron._v_scale(system)/0.673
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ v rtol = 1e-5
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ p rtol = 1E-1
        @test Clapeyron.crit_pure(system)[1]/Clapeyron.T_scale(system) ≈ 1.32 rtol = 1E-4
    end
end

@testset "SPUNG methods" begin
    system = SPUNG(["ethane"])
    p = 1e5
    T = 313.15
    @testset "Bulk properties" begin

        #GERG2008(["ethane]): 0.025868956878898026
        vv = 0.025868956878898026
        @test Clapeyron.volume(system, p, T) ≈ vv rtol = 1e-6
        @test Clapeyron.volume(system, p, T;phase=:vapour) ≈ vv rtol = 1e-6

        #GERG2008(["ethane]) : 318 m/s
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 308.38846317827625 rtol = 1e-6
    end

    T_sat = 250.15
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T_sat)[1] ≈ 1.3085074415334722e6 rtol = 1E-6
        #Critical point of ethane: 305.322
        @test Clapeyron.crit_pure(system)[1] ≈ 305.37187249327553 rtol = 1E-6
    end
end
GC.gc()
@testset "lattice methods" begin
    p = 1e5
    T = 298.15
    T1 = 301.15
    system = Clapeyron.SanchezLacombe(["carbon dioxide"])
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T1) ≈ 0.02492944175392707 rtol = 1E-6
        @test Clapeyron.speed_of_sound(system, p, T1) ≈ 307.7871016597499 rtol = 1E-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 6.468653945184592e6 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1]  ≈ 304.21081254005446 rtol = 1E-6
    end
end

@testset "PeTS" begin
    system = PeTS(["methane"])
    system.params.sigma.values[1] = 1e-10
    system.params.epsilon.values[1] = 1
    system.params.epsilon.values[1] = 1

    #Values from FeOs notebook example:
    #We can reproduce FeOs values here
    crit = crit_pure(system)
    Tc,Pc,Vc  = crit
    @test Tc ≈ 1.08905 rtol = 1e-5
    @test Vc ≈ 1/513383.86 rtol = 1e-5

    #first value from saturation pressure(T = 0.64):
    psat,vl,vv = saturation_pressure(system,0.64)
    @test psat ≈ 3.027452e+04 rtol = 1e-6
    @test vl  ≈ 1/1.359958e+06 rtol = 1e-6
    @test vv  ≈ 1/5892.917088 rtol = 1e-6

    #uses the default x0_saturation_temperature initial guess
    @test saturation_temperature(system,psat)[1] ≈ 0.64  rtol = 1e-6

    #uses the default x0_psat initial guess
    @test saturation_pressure(system,0.64,IsoFugacitySaturation(;crit))[1] ≈ 3.027452e+04 rtol = 1e-6

    T_nearc = 1.084513 #The last value of their critical point is actually above ours.
    psat_nearc = 1.374330e+06
    @test saturation_pressure(system,T_nearc)[1] ≈ psat_nearc rtol = 1e-6
    @test saturation_pressure(system,T_nearc,IsoFugacitySaturation(;crit))[1] ≈ psat_nearc rtol = 1e-6
end