using Clapeyron, Test, Unitful

@testset "SAFT methods, single components" begin
    @printline
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
        @test Clapeyron.internal_energy(system, p, T) ≈ -35882.22946560716 rtol = 1E-6
        @test Clapeyron.enthalpy(system, p, T) ≈ -35876.32155687084 rtol = 1E-6
        @test Clapeyron.gibbs_free_energy(system, p, T) ≈ -18323.87754268292 rtol = 1E-6
        @test Clapeyron.helmholtz_free_energy(system, p, T) ≈ -18329.785451419295 rtol = 1E-6
        @test Clapeyron.isochoric_heat_capacity(system, p, T) ≈ 48.37961296309505 rtol = 1E-6
        @test Clapeyron.isobaric_heat_capacity(system, p, T) ≈ 66.45719988319257 rtol = 1E-6
        @test Clapeyron.isothermal_compressibility(system, p, T) ≈ 1.1521981407243432e-9 rtol = 1E-6
        @test Clapeyron.isentropic_compressibility(system, p, T) ≈ 8.387789464951438e-10 rtol = 1E-6
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 1236.4846683094133 rtol = 1E-6 #requires that the model has Mr
        @test Clapeyron.isobaric_expansivity(system, p, T) ≈ 0.0010874255138433413 rtol = 1E-6
        @test Clapeyron.joule_thomson_coefficient(system, p, T) ≈ -6.007581864883784e-7 rtol = 1E-6
        @test Clapeyron.second_virial_coefficient(system, T) ≈ -0.004919678119638886  rtol = 1E-6 #exact value calculated by using BigFloat
        @test Clapeyron.inversion_temperature(system, 1.1e8) ≈ 824.4137805298458 rtol = 1E-6
        @test Clapeyron.fugacity_coefficient(system, p, T, phase = :l)[1] ≈ 0.07865326632570452 rtol = 1E-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 7972.550405922014 rtol = 1E-6
        @test Clapeyron.saturation_temperature(system, p)[1] ≈ 351.32529505096164 rtol = 1E-6
        @test Clapeyron.saturation_temperature(system, p, 350.)[1] ≈ 351.32529505096164 rtol = 1E-6
        @test Clapeyron.enthalpy_vap(system, T) ≈ 41712.78521121877 rtol = 1E-6
        @test Clapeyron.acentric_factor(system) ≈ 0.5730309964718605 rtol = 1E-6
        @test Clapeyron.crit_pure(system)[1] ≈ 533.1324329774004 rtol = 1E-6
    end
    @printline
end

@testset "pharmaPCSAFT, single components" begin
    system = pharmaPCSAFT(["water08"])
    v1 = Clapeyron.saturation_pressure(system, 280.15)[2]
    v2 = Clapeyron.saturation_pressure(system, 278.15)[2]
    v3 = Clapeyron.saturation_pressure(system, 275.15)[2]

    @test v1 ≈ 1.8022929328333385e-5  rtol = 1E-6
    @test v2 ≈ 1.8022662044726256e-5 rtol = 1E-6
    @test v3 ≈ 1.802442451376152e-5 rtol = 1E-6
    #density maxima of water
    @test v2 < v1
    @test v2 < v3
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
    p = 1e5
    T = 273.15 + 78.24
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T,phase=:v) ≈ 0.027368884099868623 rtol = 1e-6
        #volume(SAFTgammaMie(["ethanol"]),p,T,phase =:l)  =6.120507339375205e-5
        @test Clapeyron.volume(system, p, T,phase=:l) ≈ 6.245903786961202e-5 rtol = 1e-6
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
        #SAFT VR Mie is really sensitive to the the critical point
        @test_broken Clapeyron.crit_pure(system)[1] ≈ 522.7772913470494 rtol = 1E-5
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
        @test_broken Clapeyron.crit_pure(system)[1] ≈ 524.1501435599444  rtol = 1E-5 #TODO FIX
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

@testset "SAFT methods, multi-components" begin
    @printline
    system = PCSAFT(["methanol","cyclohexane"])
    p = 1e5
    T = 313.15
    z = [0.5,0.5]
    p2 = 2e6
    T2 = 443.15
    z2 = [0.27,0.73]
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T, z) ≈ 7.779694485714412e-5 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T, z) ≈ 1087.0303138908864 rtol = 1E-6
        @test Clapeyron.activity_coefficient(system, p, T, z)[1] ≈ 1.794138454452822 rtol = 1E-6
        @test Clapeyron.fugacity_coefficient(system, p, T, z)[1] ≈ 0.5582931304564298 rtol = 1E-6
        @test Clapeyron.mixing(system, p, T, z, Clapeyron.gibbs_free_energy) ≈ -178.10797973342596 rtol = 1E-6
        @test Clapeyron.excess(system, p, T, z, Clapeyron.volume) ≈ 1.004651584989827e-6 rtol = 1E-6
        @test Clapeyron.excess(system, p, T, z, Clapeyron.entropy) ≈ -2.832281578281112 rtol = 1E-6
        @test Clapeyron.excess(system, p, T, z, Clapeyron.gibbs_free_energy) ≈ 1626.6212908893858 rtol = 1E-6
    end
    @testset "Equilibrium properties" begin
        #Those are the highest memory-intensive routines. i suspect that this is causing the
        #failures on windows 1.6. testing if adding GC pauses helps the problem
        GC.gc()
        @test Clapeyron.gibbs_solvation(system,T) ≈ -13131.087644740426 rtol = 1E-6
        GC.gc()
        @test Clapeyron.UCEP_mix(system)[1] ≈ 319.36877456397684 rtol = 1E-6
        GC.gc()
        @test Clapeyron.bubble_pressure(system,T,z)[1] ≈ 54532.249600937736 rtol = 1E-6
        GC.gc()
        @test Clapeyron.bubble_temperature(system,p2,z)[1] ≈ 435.80890506865 rtol = 1E-6
        GC.gc()
        @test Clapeyron.dew_pressure(system,T2,z)[1] ≈ 1.6555486543884084e6 rtol = 1E-6
        GC.gc()
        @test Clapeyron.dew_temperature(system,p2,z)[1] ≈ 453.0056727580934 rtol = 1E-6
        GC.gc()
        res_LLE_p = Clapeyron.LLE_pressure(system,T,z2)
        @test res_LLE_p[1] ≈ 737971.7522006684 rtol = 1E-6
        @test Clapeyron.pressure(system,res_LLE_p[2],T,z2) ≈ Clapeyron.pressure(system,res_LLE_p[3],T,res_LLE_p[end]) rtol = 1E-6
        @test Clapeyron.pressure(system,res_LLE_p[2],T,z2) ≈ res_LLE_p[1] rtol = 1E-6
        GC.gc()
        res_LLE_T = Clapeyron.LLE_temperature(system,p,z2)
        T_LLE = res_LLE_T[1]
        @test res_LLE_T[1] ≈ 312.9523684945214  rtol = 1E-6
        @test Clapeyron.pressure(system,res_LLE_T[2],T_LLE,z2) ≈ Clapeyron.pressure(system,res_LLE_T[3],T_LLE,res_LLE_T[end]) rtol = 1E-6
        @test Clapeyron.pressure(system,res_LLE_T[2],T_LLE,z2) ≈ p rtol = 1E-6
        GC.gc()
        @test Clapeyron.azeotrope_pressure(system,T2)[1] ≈ 2.4435462800998255e6 rtol = 1E-6
        GC.gc()
        @test Clapeyron.azeotrope_temperature(system,p)[1] ≈ 328.2431049077264 rtol = 1E-6
        GC.gc()
        @test Clapeyron.UCST_mix(system,T2)[1] ≈ 1.0211532467788119e9 rtol = 1E-6
        GC.gc()
        @test Clapeyron.VLLE_pressure(system, T)[1] ≈ 54504.079665621306 rtol = 1E-6
        GC.gc()
        @test Clapeyron.VLLE_temperature(system, 54504.079665621306)[1] ≈ 313.1499860368554 rtol = 1E-6
        GC.gc()
        @test Clapeyron.crit_mix(system,z)[1] ≈ 518.0004062881115 rtol = 1E-6
    end
    @printline
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
        @test Clapeyron.wilson_k_values(system,p,T) ≈ [0.13839117786853375]  rtol = 1E-6
    end
end

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
        @test Vc == system.params.Vc.values[1]
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
        @test Pc == system.params.Pc.values[1]
        @test Vc == system.params.Vc.values[1]
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
    vc = volume(system,system.params.Pc[1],system.params.Tc[1]) #vc calculated via cubic_poly
    crit = crit_pure(system) #vc calculated via pure_cubic_zc
    @test vc ≈ crit[3] rtol = 1e-4
    @test vc/system.params.Vc[1] ≈ 1.168 rtol = 1e-4 #if Zc_exp < 0.29, this should hold, by definition
end

@testset "EPPR78, single component" begin
    system = EPPR78(["carbon dioxide"])
    T = 400u"K"
    @test Clapeyron.volume(system, 3311.0u"bar", T) ≈ 3.363141761376883e-5u"m^3"
    @test Clapeyron.molar_density(system, 3363.1u"bar", T) ≈ 29810.09484964839u"mol*m^-3"
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
        srksystem  = SRK(["ethane","undecane"])
        @test Clapeyron.wilson_k_values(srksystem,p,T) ≈ [0.420849235562207, 1.6163027384311e-5] rtol = 1E-6
    end
end

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
    if hasfield(Wilson,:puremodel)
        system2 = Wilson(["water","methanol"],puremodel = com)
    else
        system2 = CompositeModel(["water","methanol"],liquid = Wilson, fluid = com)
    end
    com1 = split_model(com)[1]
    p = 1e5
    T = 298.15
    T2 = 320.15
    z = [0.5,0.5]
    z_bulk = [0.2,0.8]

    @testset "Bulk properties" begin
        @test crit_pure(com1)[1] ≈ 647.13
        @test Clapeyron.volume(system, p, T, z_bulk) ≈ 8.602344040626639e-5 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T, z_bulk) ≈ 1371.9014493149134 rtol = 1e-6
        @test Clapeyron.mixing(system, p, T, z_bulk, Clapeyron.gibbs_free_energy) ≈ -356.86007792929263 rtol = 1e-6
        @test Clapeyron.mixing(system, p, T, z_bulk, Clapeyron.enthalpy) ≈ 519.0920708672975 rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.gibbs_solvation(system, T) ≈ -24707.145697543132 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system, T, z)[1] ≈ 23758.58099358788 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system, T, z,ActivityBubblePressure(gas_fug = true,poynting = true))[1] ≈ 23839.554959977086
        @test Clapeyron.bubble_pressure(system, T, z,ActivityBubblePressure(gas_fug = true,poynting = false))[1] ≈ 23833.324475723246

        @test Clapeyron.bubble_temperature(system,23758.58099358788, z)[1] ≈ T  rtol = 1E-6

        @test Clapeyron.dew_pressure(system2, T2, z)[1] ≈ 19386.939256733036 rtol = 1E-6
        @test Clapeyron.dew_pressure(system2, T2, z,ActivityDewPressure(gas_fug = true,poynting = true))[1] ≈ 19393.924550078184 rtol = 1e-6
        @test Clapeyron.dew_pressure(system2, T2, z,ActivityDewPressure(gas_fug = true,poynting = false))[1] ≈ 19393.76058757084 rtol = 1e-6
        @test Clapeyron.dew_temperature(system2, 19386.939256733036, z)[1]  ≈ T2 rtol = 1E-6
    end
end

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
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ 1.8068623941501927e-5 rtol = 1e-6
        @test Clapeyron.volume(system, p, T_v;phase=:vapour) ≈ 0.03116877990373624 rtol = 1e-6
        @test Clapeyron.volume(system, p_c, T_c) ≈ 0.00018553711945962424 rtol = 1e-6
        @test Clapeyron.volume(system, p_c, T_c;phase=:sc) ≈ 0.00018553711945962424 rtol = 1e-6
        @test Clapeyron.speed_of_sound(system, p, T) ≈ 1496.699163371358 rtol = 1e-6
    end
    @testset "VLE properties" begin
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 3169.9293390134403 rtol = 1E-6
        #ir varies a bit, it gives 3170.301356765357
        @test_broken Clapeyron.saturation_pressure(system, T, IsoFugacitySaturation())[1] ≈ 3169.9293390134403 rtol = 1E-6
        #saturation temperature tests are noisy
        @test Clapeyron.saturation_temperature(system,3169.9293390134403)[1] ≈ 298.1499999999789 rtol = 1E-6
        tc,pc,vc =  Clapeyron.crit_pure(system)
        @test tc ≈ 647.096 rtol = 1E-5
        v2 =  volume(system,pc,tc)
        @test pressure(system,v2,tc) ≈ pc rtol = 1E-6
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
        @test Clapeyron.saturation_pressure(system, T)[1] ≈ 97424.11102152328 rtol = 1E-6
        #they vary a litte bit. i don't know why, it gives 97423.47874065055
        @test_broken Clapeyron.saturation_pressure(system, T, IsoFugacitySaturation())[1] ≈ 97424.11102152328 rtol = 1E-6
        #saturation temperature tests are noisy
        @test Clapeyron.saturation_temperature(system,97424.11102152328)[1] ≈ 230.15014586866016  rtol = 1E-6
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

@testset "Helmholtz + Activity" begin
    model = HelmAct(["water","ethanol"])
    p = 12666.0
    x1 = Clapeyron.FractionVector( 0.00350)
    @test bubble_temperature(model,p,x1)[4][1] ≈ 0.00198 rtol = 1e-2
end

@testset "SingleFluid - CoolProp" begin
    #methanol, uses assoc term
    @test saturation_pressure(SingleFluid("methanol"),300.15)[1] ≈ PropsSI("P","T",300.15,"Q",1.,"methanol") rtol = 1e-6
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