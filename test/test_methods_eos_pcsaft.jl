@testset "SAFT methods, single components" begin
    @printline
    system = PCSAFT(["ethanol"])
    p = 1e5
    T = 298.15
    T2 = 373.15
    v = 5.907908736304141e-5
    @testset "Bulk properties" begin
        @test Clapeyron.volume(system, p, T) ≈ v rtol = 1e-6
        @test Clapeyron.volume(system, p, T;phase=:v) ≈ 0.020427920501436134 rtol = 1e-6
        @test Clapeyron.volume(system, p, T;threaded=:false) ≈ v rtol = 1e-6
        @test Clapeyron.pip(system, v, T) ≈ 6.857076349623449 rtol = 1e-6
        @test Clapeyron.is_liquid(Clapeyron.VT_identify_phase(system, v, T))
        @test Clapeyron.compressibility_factor(system, p, T) ≈ 0.002383223535444557 rtol = 1e-6
        @test Clapeyron.pressure(system, v, T) ≈ p rtol = 1e-6
        @test Clapeyron.pressure(system, 2*v, T, Clapeyron.SA[2.0]) ≈ p rtol = 1e-6
        s = Clapeyron.entropy(system, p, T)
        @test s ≈ -58.87118569239617 rtol = 1E-6
        @test Clapeyron.VT_entropy_res(system,v,T) + Clapeyron.VT_entropy(Clapeyron.idealmodel(system),v,T) ≈ s
        @test Clapeyron.chemical_potential(system, p, T)[1] ≈ -18323.877542682934 rtol = 1E-6
        u = Clapeyron.internal_energy(system, p, T)
        @test u ≈ -35882.22946560716 rtol = 1E-6
        @test Clapeyron.VT_internal_energy_res(system,v,T) + Clapeyron.VT_internal_energy(Clapeyron.idealmodel(system),v,T) ≈ u
        h = Clapeyron.enthalpy(system, p, T)
        @test h ≈ -35876.32155687084 rtol = 1E-6
        @test Clapeyron.VT_enthalpy_res(system,v,T) + Clapeyron.VT_enthalpy(Clapeyron.idealmodel(system),v,T) ≈ h
        g = Clapeyron.gibbs_free_energy(system, p, T)
        @test g ≈ -18323.87754268292 rtol = 1E-6
        @test Clapeyron.VT_gibbs_free_energy_res(system,v,T) + Clapeyron.VT_gibbs_free_energy(Clapeyron.idealmodel(system),v,T) ≈ g
        a = Clapeyron.helmholtz_free_energy(system, p, T)
        @test a ≈ -18329.785451419295 rtol = 1E-6
        @test Clapeyron.VT_helmholtz_free_energy_res(system,v,T) + Clapeyron.VT_helmholtz_free_energy(Clapeyron.idealmodel(system),v,T) ≈ a
        @test Clapeyron.isochoric_heat_capacity(system, p, T) ≈ 48.37961296309505 rtol = 1E-6
        @test Clapeyron.isobaric_heat_capacity(system, p, T) ≈ 66.45719988319257 rtol = 1E-6
        Cp = Clapeyron.isobaric_heat_capacity(system, p, T2)
        Cv = Clapeyron.isochoric_heat_capacity(system, p, T2)
        @test Clapeyron.adiabatic_index(system, p, T2) ≈ Cp/Cv rtol = 1E-12
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