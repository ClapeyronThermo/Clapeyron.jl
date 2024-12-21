@testset "Tp flash algorithms" begin
    #this is a VLLE equilibria
    system = PCSAFT(["water","cyclohexane","propane"])
    T = 298.15
    p = 1e5
    z = [0.333, 0.333, 0.334]

    @testset "RR Algorithm" begin
        method = RRTPFlash()
        @test Clapeyron.tp_flash(system, p, T, z, method)[3] ≈ -6.539976318817461 rtol = 1e-6

        #test for initialization when K suggests single phase but it could be solved supposing bubble or dew conditions.
        substances = ["water", "methanol", "propyleneglycol","methyloxirane"]
        pcp_system = PCPSAFT(substances)
        res = Clapeyron.tp_flash2(pcp_system, 25_000.0, 300.15, [1.0, 1.0, 1.0, 1.0], RRTPFlash())
        @test res.data.g ≈ -8.900576759774916 rtol = 1e-6
    end

    if isdefined(Base,:get_extension)
        @testset "RR Algorithm - MultiComponentFlash.jl" begin
            mcf = MCFlashJL()
            @test Clapeyron.tp_flash(system, p, T, z, mcf)[3] ≈ -6.490030777308265 rtol = 1e-6
        end
    end
    GC.gc()

    @testset "DE Algorithm" begin
        #VLLE eq
        @test Clapeyron.tp_flash(system, p, T, z, DETPFlash(numphases = 3))[3] ≈ -6.759674475174073 rtol = 1e-6
        #LLE eq with activities
        act_system = UNIFAC(["water","cyclohexane","propane"])
        flash0 = Clapeyron.tp_flash(act_system, p, T, [0.5,0.5,0.0], DETPFlash(equilibrium = :lle))
        act_x0 = activity_coefficient(act_system, p, T, flash0[1][1,:]) .* flash0[1][1,:]
        act_y0 = activity_coefficient(act_system, p, T, flash0[1][2,:]) .* flash0[1][2,:]
        @test Clapeyron.dnorm(act_x0,act_y0) < 0.01 #not the most accurate, but it is global
    end

    @testset "Multiphase algorithm" begin
        @test Clapeyron.tp_flash(system, p, T, z, MultiPhaseTPFlash())[3] ≈ -6.759674475175065 rtol = 1e-6
        system2 = PR(["IsoButane", "n-Butane", "n-Pentane", "n-Hexane"])
        @test Clapeyron.tp_flash(system2, 1e5, 284.4, [1,1,1,1]*0.25, MultiPhaseTPFlash())[3] ≈ -6.618441125949686 rtol = 1e-6
    end

    GC.gc()

    @testset "Michelsen Algorithm" begin
        x0 = [0.9997755902156433, 0.0002244097843566859, 0.0]
        y0 = [6.425238373915699e-6, 0.9999935747616262, 0.0]
        method = MichelsenTPFlash(x0 = x0, y0 = y0, equilibrium= :lle)
        @test Clapeyron.tp_flash(system, p, T, [0.5,0.5,0.0],method)[3] ≈ -7.577270350886795 rtol = 1e-6

        method2 = MichelsenTPFlash(x0 = x0, y0 = y0, equilibrium = :lle, ss_iters = 4, second_order = false)
        @test Clapeyron.tp_flash(system, p, T, [0.5,0.5,0.0],method2)[3] ≈ -7.577270350886795 rtol = 1e-6

        method3 = MichelsenTPFlash(x0 = x0, y0 = y0, equilibrium = :lle, ss_iters = 4,second_order = true)
        @test Clapeyron.tp_flash(system, p, T, [0.5,0.5,0.0],method3)[3] ≈ -7.577270350886795 rtol = 1e-6
    end
    GC.gc()

    @testset "Michelsen Algorithm, nonvolatiles/noncondensables support" begin

        system = PCSAFT(["hexane", "ethanol", "methane", "decane"])
        T = 320.  # K
        p = 1e5  # Pa
        z = [0.25, 0.25, 0.25, 0.25]
        x0 = [0.3, 0.3, 0., 0.4]
        y0 = [0.2, 0.2, 0.6, 0.]

        method_normal = MichelsenTPFlash(x0=x0, y0=y0, second_order=true)
        @test Clapeyron.tp_flash(system, p, T, z, method_normal)[1] ≈
        [0.291924  0.306002  0.00222251  0.399851
        0.181195  0.158091  0.656644    0.00406898] rtol = 1e-6
        GC.gc()

        method_nonvolatiles = MichelsenTPFlash(x0=x0, y0=y0, ss_iters = 1, second_order=true, nonvolatiles = ["decane"])
        @test Clapeyron.tp_flash(system, p, T, z, method_nonvolatiles)[1] ≈
        [0.291667  0.305432  0.00223826  0.400663
        0.180861  0.15802   0.661119    0.0] rtol = 1e-6
        GC.gc()

        method_noncondensables = MichelsenTPFlash(x0=x0, y0=y0,ss_iters = 1, second_order=false, noncondensables = ["methane"])
        @test Clapeyron.tp_flash(system, p, T, z, method_noncondensables)[1] ≈
        [0.292185  0.306475  0.0      0.40134
        0.181452  0.158233  0.65623  0.00408481] rtol = 1e-6
        GC.gc()

        method_both = MichelsenTPFlash(x0=x0, y0=y0, ss_iters = 1,second_order=false, noncondensables = ["methane"],nonvolatiles = ["decane"])
        @test Clapeyron.tp_flash(system, p, T, z, method_both)[1] ≈
        [0.291928  0.3059    0.0       0.402171
        0.181116  0.158162  0.660722  0.0] rtol = 1e-6
        GC.gc()
    end

    @testset "Michelsen Algorithm, activities" begin
    #example from https://github.com/ClapeyronThermo/Clapeyron.jl/issues/144
        system = UNIFAC(["water", "hexane"])
        alg1 = MichelsenTPFlash(
            equilibrium = :lle,
            K0 = [0.00001/0.99999, 0.99999/0.00001],
        )

        flash1 = tp_flash(system, 101325, 303.15, [0.5, 0.5], alg1)
        act_x1 = activity_coefficient(system, 101325, 303.15, flash1[1][1,:]) .* flash1[1][1,:]
        act_y1 = activity_coefficient(system, 101325, 303.15, flash1[1][2,:]) .* flash1[1][2,:]
        @test Clapeyron.dnorm(act_x1,act_y1) < 1e-8

        alg2 = RRTPFlash(
            equilibrium = :lle,
            x0 = [0.99999, 0.00001],
            y0 = [0.00001, 0.00009]
        )
        flash2 = tp_flash(system, 101325, 303.15, [0.5, 0.5], alg2)
        act_x2 = activity_coefficient(system, 101325, 303.15, flash2[1][1,:]) .* flash2[1][1,:]
        act_y2 = activity_coefficient(system, 101325, 303.15, flash2[1][2,:]) .* flash2[1][2,:]
        @test Clapeyron.dnorm(act_x2,act_y2) < 1e-8

        #test K0_lle_init initialization
        alg3 = RRTPFlash(
            equilibrium = :lle)
        flash3 = tp_flash(system, 101325, 303.15, [0.5, 0.5], alg3)
        act_x3 = activity_coefficient(system, 101325, 303.15, flash3[1][1,:]) .* flash3[1][1,:]
        act_y3 = activity_coefficient(system, 101325, 303.15, flash3[1][2,:]) .* flash3[1][2,:]
        @test Clapeyron.dnorm(act_x3,act_y3) < 1e-8

        #test combinations of Activity + CompositeModel
        system_fluid = CompositeModel(["water","ethanol"],gas = BasicIdeal, liquid = RackettLiquid, saturation = LeeKeslerSat)
        system_cc  = CompositeModel(["water","ethanol"],liquid = UNIFAC,fluid = system_fluid)
        flash3 = tp_flash(system_cc, 101325, 303.15, [0.5, 0.5], alg2)
        act_x3 = activity_coefficient(system_cc, 101325, 303.15, flash3[1][1,:]) .* flash3[1][1,:]
        act_y3 = activity_coefficient(system_cc, 101325, 303.15, flash3[1][2,:]) .* flash3[1][2,:]
        @test Clapeyron.dnorm(act_x3,act_y3) < 1e-8

        #running the vle part
        if hasfield(UNIFAC,:puremodel)
            model_vle = UNIFAC(["water", "ethanol"],puremodel = PCSAFT)
        else
            model_vle = CompositeModel(["water", "ethanol"],liquid = UNIFAC,fluid = PCSAFT)
        end
        flash4 = tp_flash(model_vle, 101325, 363.15, [0.5, 0.5], MichelsenTPFlash())
        #=@test flash4[1] ≈
        [0.6824441505154921 0.31755584948450793
        0.3025308123759482 0.6974691876240517] rtol = 1e-6
        this was wrong, we were calculating the gas volume as the addition of partial pressures,
        basically ideal gas.
        =#

        @test flash4[1] ≈
        [0.7006206854062672 0.29937931459373285;
        0.43355504959745633 0.5664449504025437] rtol = 1e-6
        #test equality of activities does not make sense in VLE
    end

    @testset "Michelsen Algorithm, CompositeModel" begin
        p,T,z = 101325.,85+273.,[0.2,0.8]
        system = CompositeModel(["water","ethanol"],gas = BasicIdeal, liquid = RackettLiquid, saturation = LeeKeslerSat) #ideal gas + rackett + lee kesler saturation correlation
        @test Clapeyron.tp_flash(system, p, T, z, MichelsenTPFlash())[1] ≈
        [0.3618699659002134 0.6381300340997866
        0.17888243361092543 0.8211175663890746] rtol = 1e-6

        @test_throws ErrorException Clapeyron.tp_flash(system, p, T, z, MichelsenTPFlash(ss_iters = 0))
    end
end

@testset "XY flash" begin
    #1 phase (#320)
    model = cPR(["ethane","methane"],idealmodel = ReidIdeal)
    p = 101325.0
    z = [1.2,1.2]
    T = 350.0
    h = enthalpy(model,p,T,z)
    res0 = ph_flash(model,p,h,z)
    @test Clapeyron.temperature(res0) ≈ T rtol = 1e-6
    @test enthalpy(model,res0) ≈ h rtol = 1e-6

    #2 phases
    h = -13831.0
    res1 = ph_flash(model,p,h,z)
    @test enthalpy(model,res1) ≈ h rtol = 1e-6

    #test for ps_flash:
    model = cPR(["ethane"],idealmodel = ReidIdeal);
    s = 100;p=101325;z = [1.0];
    h = Clapeyron.PS.enthalpy(model,p,s,z)
    s2 = Clapeyron.PH.entropy(model,p,h,z)
    @test s ≈ s2 rtol = 1e-6
    #examples for qt, qp flash (#314)
    model = cPR(["ethane","propane"],idealmodel=ReidIdeal)
    res2 = qt_flash(model,0.5,208.0,[0.5,0.5])
    @test Clapeyron.pressure(res2) ≈ 101634.82435966855 rtol = 1e-6
    res3 = qp_flash(model,0.5,120000.0,[0.5,0.5])
    @test Clapeyron.temperature(res3) ≈ 211.4972567716822 rtol = 1e-6

    #1 phase input should error
    model = PR(["IsoButane", "n-Butane", "n-Pentane", "n-Hexane"])
    z = [0.25, 0.25, 0.25, 0.25]
    p = 1e5
    h = enthalpy(model, 1e5, 303.15, z)
    r = Clapeyron.ph_flash(model, p, h, z)
    @test_throws ArgumentError qt_flash(model,0.5,308,z,flash_result = r)
    res4 = qp_flash(model,0.7,60000.0,z)
    T4 = Clapeyron.temperature(res4)
    @test pressure(model,res4.volumes[1],T4,res4.compositions[1]) ≈ 60000.0 rtol = 1e-6
    @test pressure(model,res4.volumes[2],T4,res4.compositions[2]) ≈ 60000.0 rtol = 1e-6

    #example in documentation for xy_flash
    spec = FlashSpecifications(p = 101325.0, T = 200.15) #p-T flash
    model = cPR(["ethane","propane"],idealmodel=ReidIdeal)
    z = [0.5,0.5] #bulk composition
    x1 = [0.25,0.75] #liquid composition
    x2 = [0.75,0.25] #gas composition
    compositions = [x1,x2]
    volumes = [6.44e-5,0.016]
    fractions = [0.5,0.5]
    p0,T0 = NaN,NaN #in p-T flash, pressure and temperature are already specifications
    data = FlashData(p0,T0)
    result0 = FlashResult(compositions,fractions,volumes,data) #a FlashResult containing all necessary information
    result = xy_flash(model,spec,z,result0) #perform the flash
    @test Clapeyron.temperature(result) == 200.15
    @test pressure(result) == 101325.0

    #px_flash_pure/tx_flash_pure, 1 phase (#320)
    model = cPR(["ethane"],idealmodel = ReidIdeal);
    p = 101325;h = 100; z = Clapeyron.SA[1]; T = Clapeyron.PH.temperature(model,p,h,z)
    @test enthalpy(model,p,T,z) ≈ h rtol = 1e-6
    res5 = Clapeyron.tx_flash_pure(model,T,h,z,enthalpy)
    @test pressure(res5)  ≈ p rtol = 1e-6
end

@testset "Saturation Methods" begin
    model = PR(["water"])
    vdw = vdW(["water"])
    p0 = 1e5
    T = 373.15
    p,vl,vv = Clapeyron.saturation_pressure(model,T) #default

    #legacy api,
    @test Clapeyron.saturation_pressure(model,T,Clapeyron.ChemPotVSaturation((vl,vv)))[1] ==
        Clapeyron.saturation_pressure(model,T,Clapeyron.ChemPotVSaturation([vl,vv]))[1] ==
        Clapeyron.saturation_pressure(model,T,[vl,vv])[1] ==
        Clapeyron.saturation_pressure(model,T,(vl,vv))[1]

    px,vlx,vvx = Clapeyron.saturation_pressure(vdw,T) #vdw

    p1,vl1,vv1 = Clapeyron.saturation_pressure_impl(model,T,IsoFugacitySaturation())
    @test p1 ≈ p rtol = 1e-6
    p2,vl2,vv2 = Clapeyron.saturation_pressure_impl(model,T,IsoFugacitySaturation(p0 = 1e5))
    @test p1 ≈ p rtol = 1e-6
    p3,vl3,vv3 = Clapeyron.saturation_pressure_impl(model,T,ChemPotDensitySaturation())
    @test p3 ≈ p rtol = 1e-6
    p4,vl4,vv4 = Clapeyron.saturation_pressure_impl(model,T,ChemPotDensitySaturation(;vl,vv))
    p4b,vl4b,vv4b = Clapeyron.psat_chempot(model,T,vl,vv)
    @test p4 ≈ p rtol = 1e-6
    @test (p4 == p4b) && (vl4 == vl4b) && (vv4 == vv4b)
    GC.gc()

    #test IsoFugacity, near criticality
    Tc_near = 0.95*647.096
    psat_Tcnear = 1.496059652088857e7 #default solver result
    if Base.VERSION < v"1.11" && Sys.islinux()
        @test first(Clapeyron.saturation_pressure(model,Tc_near,IsoFugacitySaturation(p0 = 1.49e7))) ≈ psat_Tcnear rtol = 1e-6
    else
        @test first(Clapeyron.saturation_pressure(model,Tc_near,IsoFugacitySaturation())) ≈ psat_Tcnear rtol = 1e-6

    end
    #Test that IsoFugacity fails over critical point
    @test isnan(first(Clapeyron.saturation_pressure(model,1.1*647.096,IsoFugacitySaturation())))
    GC.gc()

    #SuperAncSaturation
    p5,vl5,vv5 = Clapeyron.saturation_pressure_impl(model,T,SuperAncSaturation())
    @test p5 ≈ p rtol = 1e-6
    @test Clapeyron.saturation_temperature_impl(model,p5,SuperAncSaturation())[1] ≈ T rtol = 1e-6
    @test @inferred Clapeyron.saturation_pressure_impl(vdw,T,SuperAncSaturation())[1] ≈ px
    GC.gc()

    #AntoineSat
    @test Clapeyron.saturation_temperature(model,p0,AntoineSaturation(T0 = 400.0))[1] ≈ 374.2401401001685 rtol = 1e-6
    @test Clapeyron.saturation_temperature(model,p0,AntoineSaturation(vl = vl5,vv = vv5))[1] ≈ 374.2401401001685 rtol = 1e-6
    @test_throws Any Clapeyron.saturation_temperature(model,p0,AntoineSaturation(vl = vl5,T0 = 400))
    GC.gc()

    #ClapeyronSat
    @test Clapeyron.saturation_temperature(model,p0,ClapeyronSaturation())[1] ≈ 374.2401401001685 rtol = 1e-6

    #Issue #290
    @test Clapeyron.saturation_temperature(cPR("R1233zde"),101325*20,crit_retry = false)[1] ≈ 405.98925205830335 rtol = 1e-6
    @test Clapeyron.saturation_temperature(cPR("isobutane"),1.7855513185537157e6,crit_retry = false)[1] ≈ 366.52386488214876 rtol = 1e-6
    @test Clapeyron.saturation_temperature(cPR("propane"),2.1298218093361156e6,crit_retry = false)[1] ≈ 332.6046103831853 rtol = 1e-6
    @test Clapeyron.saturation_temperature(cPR("r134a"),2.201981727901889e6,crit_retry = false)[1] ≈ 344.6869001549851 rtol = 1e-6
end

@testset "Tproperty" begin
    model1 = PCSAFT(["propane","dodecane"])
    p = 101325.0; T = 300.0;z = [0.5,0.5]
    h_ = enthalpy(model1,p,T,z)
    s_ = entropy(model1,p,T,z)
    @test Tproperty(model1,p,h_,z,enthalpy) ≈ T
    @test Tproperty(model1,p,s_,z,entropy) ≈ T

    model2 = PCSAFT(["propane"])
    z2 = [1.]
    h2_ = enthalpy(model2,p,T,z2)
    s2_ = entropy(model2,p,T,z2)
    @test Tproperty(model2,p,h2_,z2,enthalpy) ≈ T
    @test Tproperty(model2,p,s2_,z2,entropy) ≈ T

    #issue 309
    model3 = cPR(["ethane"],idealmodel=ReidIdeal)
    T3 = 300
    z3 = [5]
    s30 = entropy(model3,p,T3,z3)
    p3 = 2*p
    T3_calc = Tproperty(model3,p3,s30,z3,entropy)
    s3 = entropy(model3,p3,T3_calc,z3)
    @test s3 ≈ s30

    #issue 309 (https://github.com/ClapeyronThermo/Clapeyron.jl/issues/309#issuecomment-2508038968)
    model4 = cPR("R134A",idealmodel= ReidIdeal)
    T_crit,p_crit,_ = crit_pure(model4)
    T1 = 300.0
    p1 = saturation_pressure(model4,T1)[1] + 101325
    s1 = entropy(model4,p1,T1)
    h1 = enthalpy(model4,p1,T1)
    p2 = p_crit + 2*101325
    T2 =  Tproperty(model4,p2,s1,Clapeyron.SA[1.0],entropy)
    s2 = entropy(model4,p2,T2)
    h2 = enthalpy(model4,p2,T2)
    @test s2 ≈ s1
end

@testset "bubble/dew point algorithms" begin
    system1 = PCSAFT(["methanol","cyclohexane"])
    p = 1e5
    T = 313.15
    z = [0.5,0.5]
    p2 = 2e6
    T2 = 443.15
    z2 = [0.27,0.73]

    pres1 = 54532.249600937736
    Tres1 = 435.80890506865
    pres2 = 1.6555486543884084e6
    Tres2 = 453.0056727580934
    @testset "bubble pressure" begin
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.ChemPotBubblePressure())[1] ≈ pres1 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.ChemPotBubblePressure(y0 = [0.6,0.4]))[1] ≈ pres1 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.ChemPotBubblePressure(p0 = 5e4))[1] ≈ pres1 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.ChemPotBubblePressure(p0 = 5e4,y0 = [0.6,0.4]))[1] ≈ pres1 rtol = 1E-6
        GC.gc()

        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.FugBubblePressure())[1] ≈ pres1 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.FugBubblePressure(y0 = [0.6,0.4]))[1] ≈ pres1 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.FugBubblePressure(p0 = 5e4))[1] ≈ pres1 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.FugBubblePressure(p0 = 5e4,y0 = [0.6,0.4]))[1] ≈ pres1 rtol = 1E-6
        #test multidimensional fugacity solver
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.FugBubblePressure(itmax_newton = 1))[1]  ≈ pres1 rtol = 1E-6
        GC.gc()

        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.ActivityBubblePressure())[1] ≈ pres1 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.ActivityBubblePressure(y0 = [0.6,0.4]))[1] ≈ pres1 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.ActivityBubblePressure(p0 = 5e4))[1] ≈ pres1 rtol = 1E-6
        @test Clapeyron.bubble_pressure(system1,T,z,Clapeyron.ActivityBubblePressure(p0 = 5e4,y0 = [0.6,0.4]))[1] ≈ pres1 rtol = 1E-6
        GC.gc()
    end

    @testset "bubble temperature" begin
        @test Clapeyron.bubble_temperature(system1,p2,z,Clapeyron.ChemPotBubbleTemperature())[1] ≈ Tres1 rtol = 1E-6
        @test Clapeyron.bubble_temperature(system1,p2,z,Clapeyron.ChemPotBubbleTemperature(y0 = [0.7,0.3]))[1] ≈ Tres1 rtol = 1E-6
        @test Clapeyron.bubble_temperature(system1,p2,z,Clapeyron.ChemPotBubbleTemperature(T0 = 450))[1] ≈ Tres1 rtol = 1E-6
        @test Clapeyron.bubble_temperature(system1,p2,z,Clapeyron.ChemPotBubbleTemperature(T0 = 450,y0 = [0.75,0.25]))[1] ≈ Tres1 rtol = 1E-6
        GC.gc()

        @test Clapeyron.bubble_temperature(system1,p2,z,Clapeyron.FugBubbleTemperature())[1] ≈ Tres1 rtol = 1E-6
        @test Clapeyron.bubble_temperature(system1,p2,z,Clapeyron.FugBubbleTemperature(y0 = [0.75,0.25]))[1] ≈ Tres1 rtol = 1E-6
        @test Clapeyron.bubble_temperature(system1,p2,z,Clapeyron.FugBubbleTemperature(T0 = 450))[1] ≈ Tres1 rtol = 1E-6
        @test Clapeyron.bubble_temperature(system1,p2,z,Clapeyron.FugBubbleTemperature(T0 = 450,y0 = [0.75,0.25]))[1] ≈ Tres1 rtol = 1E-6
        @test Clapeyron.bubble_temperature(system1,p2,z,Clapeyron.FugBubbleTemperature(itmax_newton = 1))[1] ≈ Tres1 rtol = 1E-6
        GC.gc()
    end

    @testset "dew pressure" begin
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.ChemPotDewPressure())[1] ≈ pres2 rtol = 1E-6
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.ChemPotDewPressure(x0 = [0.1,0.9]))[1] ≈ pres2 rtol = 1E-6
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.ChemPotDewPressure(p0 = 1.5e6))[1] ≈ pres2 rtol = 1E-6
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.ChemPotDewPressure(p0 = 1.5e6,x0 = [0.1,0.9]))[1] ≈ pres2 rtol = 1E-6
        GC.gc()

        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.FugDewPressure())[1] ≈ pres2 rtol = 1E-6
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.FugDewPressure(x0 = [0.1,0.9]))[1] ≈ pres2 rtol = 1E-6
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.FugDewPressure(p0 = 1.5e6))[1] ≈ pres2 rtol = 1E-6
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.FugDewPressure(p0 = 1.5e6,x0 = [0.1,0.9]))[1] ≈ pres2 rtol = 1E-6
        #for some reason, it requires 2 newton iterations.
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.FugDewPressure(itmax_newton = 2))[1] ≈ pres2 rtol = 1E-6
        GC.gc()
        #not exactly the same results, as activity coefficients are ultimately an aproximation of the real helmholtz function.
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.ActivityDewPressure())[1] ≈ pres2 rtol = 1E-3
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.ActivityDewPressure(x0 = [0.1,0.9]))[1] ≈ pres2 rtol = 1E-3
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.ActivityDewPressure(p0 = 1.5e6))[1] ≈ pres2 rtol = 1E-3
        @test Clapeyron.dew_pressure(system1,T2,z,Clapeyron.ActivityDewPressure(p0 = 1.5e6,x0 = [0.1,0.9]))[1] ≈ pres2 rtol = 1E-3
        GC.gc()
    end

    @testset "dew temperature" begin
        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.ChemPotDewTemperature())[1] ≈ Tres2 rtol = 1E-6
        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.ChemPotDewTemperature(x0 = [0.1,0.9]))[1] ≈ Tres2 rtol = 1E-6
        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.ChemPotDewTemperature(T0 = 450))[1] ≈ Tres2 rtol = 1E-6
        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.ChemPotDewTemperature(T0 = 450,x0 = [0.1,0.9]))[1] ≈ Tres2 rtol = 1E-6
        GC.gc()

        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.FugDewTemperature())[1] ≈ Tres2 rtol = 1E-6
        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.FugDewTemperature(x0 = [0.1,0.9]))[1] ≈ Tres2 rtol = 1E-6
        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.FugDewTemperature(T0 = 450))[1] ≈ Tres2 rtol = 1E-6
        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.FugDewTemperature(T0 = 450,x0 = [0.1,0.9]))[1] ≈ Tres2 rtol = 1E-6
        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.FugDewTemperature(itmax_newton = 2))[1] ≈ Tres2 rtol = 1E-6
        GC.gc()
    end

    #nonvolatiles/noncondensables testing. it also test model splitting
    system2 = PCSAFT(["hexane", "ethanol", "methane", "decane"])
    T = 320.  # K
    p = 1e5  # Pa
    z = [0.25, 0.25, 0.25, 0.25]
    x0 = [0.3, 0.3, 0., 0.4]
    y0 = [0.2, 0.2, 0.6, 0.]
    pres1 = 33653.25605767739
    Tres1 = 349.3673410368543
    pres2 = 112209.1535730352
    Tres2 = 317.58287413031866

    @testset "bubble pressure - nonvolatiles" begin
        (pa,vla,vva,ya) = bubble_pressure(system2,T,x0,FugBubblePressure(y0 = y0,p0 = 1e5,nonvolatiles = ["decane"]))
        @test pa  ≈ pres1 rtol = 1E-6
        @test ya[4] == 0.0
        (pb,vlb,vvb,yb) = bubble_pressure(system2,T,x0,ChemPotBubblePressure(y0 = y0,p0 = 1e5,nonvolatiles = ["decane"]))
        @test pa  ≈ pres1 rtol = 1E-6
        @test ya[4] == 0.0
    end
    GC.gc()

    @testset "bubble temperature - nonvolatiles" begin
        (Ta,vla,vva,ya) = bubble_temperature(system2,p,x0,FugBubbleTemperature(y0 = y0,T0 = T,nonvolatiles = ["decane"]))
        @test Ta  ≈ Tres1 rtol = 1E-6
        @test ya[4] == 0.0
        (Tb,vlb,vvb,yb) = bubble_temperature(system2,p,x0,ChemPotBubbleTemperature(y0 = y0,T0 = T,nonvolatiles = ["decane"]))
        @test Tb  ≈ Tres1 rtol = 1E-6
        @test yb[4] == 0.0
    end
    GC.gc()

    @testset "dew pressure - noncondensables" begin
        (pa,vla,vva,xa) = dew_pressure(system2,T,y0,FugDewPressure(noncondensables = ["methane"],p0 = p,x0 = x0))
        @test pa  ≈ pres2 rtol = 1E-6
        @test xa[3] == 0.0
        (pb,vlb,vvb,xb) = dew_pressure(system2,T,y0,ChemPotDewPressure(noncondensables = ["methane"],p0 = p,x0 = x0))
        @test pb  ≈ pres2 rtol = 1E-6
        @test xa[3] == 0.0
    end
    GC.gc()

    @testset "dew temperature - noncondensables" begin
        (Ta,vla,vva,xa) = dew_temperature(system2,p,y0,FugDewTemperature(noncondensables = ["methane"],T0 = T,x0 = x0))
        @test Ta  ≈ Tres2 rtol = 1E-6
        @test xa[3] == 0.0
        (Tb,vlb,vvb,xb) = dew_temperature(system2,p,y0,ChemPotDewTemperature(noncondensables = ["methane"],T0 = T,x0 = x0))
        @test Tb  ≈ Tres2 rtol = 1E-6
        @test xa[3] == 0.0
    end
    GC.gc()
end
