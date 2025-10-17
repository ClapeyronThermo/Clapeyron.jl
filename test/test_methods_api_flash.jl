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

        #https://julialang.zulipchat.com/#narrow/channel/265161-Clapeyron.2Ejl/topic/The.20meaning.20of.20subcooled.20liquid.20flash.20results
        z_zulip1 = [0.25, 0.25, 0.25, 0.25]
        p_zulip1 = 1e5
        model_zulip1 = PR(["IsoButane", "n-Butane", "n-Pentane", "n-Hexane"])
        #bubble_temperature(model, p, z) # 282.2827723244425 K
        res1 = Clapeyron.tp_flash2(model_zulip1, p_zulip1, 282.2, z_zulip1, RRTPFlash(equilibrium=:vle))
        res2 = Clapeyron.tp_flash2(model_zulip1, p_zulip1, 282.3, z_zulip1, RRTPFlash(equilibrium=:vle))
        if Clapeyron.numphases(res1) == 2
            @test iszero(res1.fractions[2])
            @test res1.volumes[1] ≈ 0.00010665596678830227 rtol = 1e-6
        else
            @test res1.volumes[1] ≈ 0.00010665596678830227 rtol = 1e-6
        end

        @test Clapeyron.numphases(res2) == 2
        @test res2.fractions[1] ≈ 0.9991083897702745 rtol = 1e-6

        #https://julialang.zulipchat.com/#narrow/channel/265161-Clapeyron.2Ejl/topic/The.20meaning.20of.20subcooled.20liquid.20flash.20results/near/534216551
        model_zulip2 = PR(["n-butane", "n-pentane", "n-hexane", "n-heptane"])
        res3 = Clapeyron.tp_flash2(model_zulip2, 1e5 , 450, z_zulip1, RRTPFlash(equilibrium=:vle))
        
        if Clapeyron.numphases(res3) == 2
            @test isone(res3.fractions[2])
            @test res3.volumes[1] ≈ 0.03683358805181434 rtol = 1e-6
        else
            @test res3.volumes[1] ≈ 0.03683358805181434 rtol = 1e-6
        end
    end

    if isdefined(Base,:get_extension)
        @testset "MultiComponentFlash.jl Algorithm" begin

            #two-phase test, using Clapeyron api
            mcf = MCFlashJL()
            @test Clapeyron.numphases(Clapeyron.tp_flash2(system, p, T, z, mcf)) == 2
            #vapour test, using MCF api
            cond = (p = 5e6, T = 303.15, z = [0.4, 0.6])
            vapour_model = PR78(["hydrogen", "methane"])
            vapour_res = MultiComponentFlash.flashed_mixture_2ph(vapour_model,cond)
            @test vapour_res.state == MultiComponentFlash.single_phase_v
            @test vapour_res.vapor.Z ≈ 0.9672507136048648 rtol = 1e-6

            #liquid test,using MCF api
            liquid_model = cPR(["octane","nonane"])
            cond = (p = 5e7, T = 303.15, z = [0.4, 0.6])
            liquid_res = MultiComponentFlash.flashed_mixture_2ph(liquid_model,cond)
            @test liquid_res.state == MultiComponentFlash.single_phase_l
            @test liquid_res.liquid.Z ≈ 3.458550315299117 rtol = 1e-6
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
        res0 = Clapeyron.tp_flash2(system, p, T, [0.5,0.5,0.0],method)
        @test Clapeyron.tp_flash(system, p, T, [0.5,0.5,0.0],method)[3] ≈ -7.577270350886795 rtol = 1e-6
        @test Clapeyron.tp_flash(system,p,T,[0.5,0.5,0.0], MichelsenTPFlash(flash_result = res0,equilibrium = :lle))[3] ≈ -7.577270350886795 rtol = 1e-6
        method2 = MichelsenTPFlash(x0 = x0, y0 = y0, equilibrium = :lle, ss_iters = 4, second_order = false)
        @test Clapeyron.tp_flash(system, p, T, [0.5,0.5,0.0],method2)[3] ≈ -7.577270350886795 rtol = 1e-6

        method3 = MichelsenTPFlash(x0 = x0, y0 = y0, equilibrium = :lle, ss_iters = 4,second_order = true)
        @test Clapeyron.tp_flash(system, p, T, [0.5,0.5,0.0],method3)[3] ≈ -7.577270350886795 rtol = 1e-6

        @testset "#454" begin
            mix = PR(["n-butane", "n-pentane", "n-hexane", "n-heptane"];
                        idealmodel=AlyLeeIdeal,
                        userlocations=(;
                            Tc             = [425.12, 469.7, 507.6, 540.2],
                            Pc             = [37.96e5, 33.7e5, 30.25e5, 27.4e5],
                            Mw             = [58.1234, 72.15028, 86.17716, 100.20404],
                            acentricfactor = [0.200164, 0.251506, 0.301261, 0.349469],
                            k              = [
                            0.0        0.0174       -0.0056      0.0033
                            0.0174     0.0          -0.00071726  0.0074
                            -0.0056    -0.00071726    0.0        -0.0078
                            0.0033     0.0074       -0.0078      0.0],
                            l              = zeros(4, 4)
                        )
                    )

            res1 = Clapeyron.tp_flash2(mix, 153_823.0, 321.9670623578307, [0.007682, 0.9923, 1.517e-17, 1.918e-31], RRTPFlash(equilibrium = :vle))
            @test res1.compositions[1] ≈ [0.0023666624484214222, 0.9976333375515787, 0.0, 0.0] rtol = 1e-6

            res2 = Clapeyron.tp_flash2(mix, 701739.83, 430.74, [2.984e-14, 0.0615, 3.48, 2.059], RRTPFlash(equilibrium = :vle))
            @test res2.compositions[1] ≈ [5.306960867808201e-15, 0.010962897986743346, 0.6212133688148559, 0.3678237331983954] rtol = 1e-6

            res3 = Clapeyron.tp_flash2(mix, 1.985550610608908e6, 416.6628781711617, [55.461373286206445, 0.09264900343401582, 7.265116936961075e-9, 8.855321114218425e-14], RRTPFlash(equilibrium = :vle))
            @test iszero(res3.fractions[1])
            @test res3.fractions[2] ≈ 55.554022296905664

            res4 = Clapeyron.tp_flash2(mix, 5.35202e5, 393.265, [36.495044786426966, 0.005798955283355085, 1.9416516061189107e-10, 2.0015179988524742e-15], RRTPFlash(equilibrium=:vle, verbose=true))
            @test iszero(res4.fractions[1])
            @test res4.fractions[2] ≈ 36.50084374190448

            res5 = Clapeyron.tp_flash2(mix, 442595.31887270656, 318.91991913774194, [18.697907101753938, 9.208988950434023e-8, 2.317361697667793e-22, 1.9317538045050555e-32], RRTPFlash(equilibrium=:vle, verbose=true))
            @test res5.fractions[1] ≈ 18.69790719384074 rtol = 1e-4
        end

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

        #water-oxygen system, non-condensables
        model_a_ideal = CompositeModel(["water","oxygen"],liquid = RackettLiquid,gas = BasicIdeal,saturation = DIPPR101Sat)
        @test Clapeyron.tp_flash(model_a_ideal,134094.74892634258,70 + 273.15,[18500.0, 24.08],noncondensables = ["oxygen"])[1] ≈
        [1.0 0.0;
        0.23252954843762222 0.7674704515623778] rtol = 1e-6

        #403
        model403 = PCSAFT(["water","carbon dioxide"])
        res = Clapeyron.tp_flash2(model403, 1e5, 323.15,[0.5,0.5],MichelsenTPFlash(nonvolatiles=["water"]))
        @test res.compositions[2] == [0.,1.]
        @test res.compositions[1] ≈ [0.999642, 0.000358065] rtol = 1e-6
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
        system_fluid = CompositeModel(["water", "hexane"],gas = BasicIdeal, liquid = RackettLiquid, saturation = LeeKeslerSat)
        system_cc  = CompositeModel(["water", "hexane"],liquid = UNIFAC,fluid = system_fluid)
        flash3 = tp_flash(system_cc, 101325, 303.15, [0.5, 0.5], alg2)
        act_x3 = activity_coefficient(system_cc, 101325, 303.15, flash3[1][1,:]) .* flash3[1][1,:]
        act_y3 = activity_coefficient(system_cc, 101325, 303.15, flash3[1][2,:]) .* flash3[1][2,:]
        @test Clapeyron.dnorm(act_x3,act_y3) < 1e-8

        #running the vle part
        if hasfield(UNIFAC,:puremodel)
            model_vle = UNIFAC(["octane","heptane"],puremodel = cPR)
        else
            model_vle = CompositeModel(["octane","heptane"],liquid = UNIFAC,fluid = cPR)
        end
        flash4 = tp_flash(model_vle, 2500.0, 300.15, [0.9, 0.1], MichelsenTPFlash())

        @test flash4[1] ≈
        [0.923964726801428 0.076035273198572;
        0.7934765930306608 0.20652340696933932] rtol = 1e-6
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
    @test PH.temperature(model,p,h,z) ≈ T rtol = 1e-6
    @test Clapeyron.temperature(PH.flash(model,p,h,z)) ≈ T rtol = 1e-6
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
    @test QT.pressure(model,0.5,208.0,[0.5,0.5]) ≈ 101634.82435966855 rtol = 1e-6
    res3 = qp_flash(model,0.5,120000.0,[0.5,0.5])
    @test Clapeyron.temperature(res3) ≈ 211.4972567716822 rtol = 1e-6
    @test QP.temperature(model,0.5,120000.0,[0.5,0.5]) ≈ 211.4972567716822 rtol = 1e-6
    #1 phase input should error
    model = PR(["IsoButane", "n-Butane", "n-Pentane", "n-Hexane"])
    z = [0.25, 0.25, 0.25, 0.25]
    p = 1e5
    h = 6300.0
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
    @test pressure(res5) ≈ p rtol = 1e-6

    #px_flash_pure: two phase (#320)
    model = cPR(["ethane"],idealmodel = ReidIdeal)
    p = 101325; z = [5.0];
    T = saturation_temperature(model,p)[1]
    h_liq = enthalpy(model,p,T-0.1,z);h_gas = enthalpy(model,p,T+0.1,z)
    h = (h_liq + h_gas)/2
    T2 = Clapeyron.PH.temperature(model,p,h,z)
    @test T2 ≈ T rtol = 1e-6

    #QP flash pure (#325)
    model = cPR(["methane"],idealmodel = ReidIdeal)
    p = 101325.0; q = 1.0; z = [5.0]
    res_qp1 = qp_flash(model,q,p,z)
    @test Clapeyron.temperature(res_qp1) == saturation_temperature(model,p)[1]
    res_qt1 = qt_flash(model,0.99,160.0,z)
    @test pressure(res_qt1) == saturation_pressure(model,160.0)[1]

    #QP: test a range of vapour fractions
    q = range(0.001,0.999,100)
    fluids = ["isobutane","pentane"]
    model = cPR(fluids,idealmodel=ReidIdeal)
    p = 101325.0
    z = [5.0,5.0]
    qp_flashes = qp_flash.(model,q,p,Ref(z))
    T = Clapeyron.temperature.(qp_flashes)
    @test maximum(diff(T)) < 0.25

    #bubble/dew temperatures via qp_flash
    Tbubble0 = bubble_temperature(model,p,z)[1]
    Tbubble1 = Clapeyron.temperature(qp_flash(model,0,p,z))
    @test Tbubble0 ≈ Tbubble1 rtol = 1e-6

    Tdew0 = dew_temperature(model,p,z)[1]
    Tdew1 = Clapeyron.temperature(qp_flash(model,1,p,z))
    #somehow, this test only fails in ubuntu-latest 1.11.3
    #@test Tdew0 ≈ Tdew1 rtol = 1e-6

    #qp_flash unmasked an error in the calculation of the initial K-values (#325)
    fluids = ["isobutane","toluene"]
    model = cPR(fluids,idealmodel=ReidIdeal)
    p = 8*101325.0; z = [5.0,5.0];
    res_qp2 = qp_flash(model,0.4,p,z)
    @test res_qp2.fractions ≈ [6.0,4.0]

    #qp_flash scaling error (#325)
    fluids = ["isopentane","isobutane"]
    model = cPR(fluids,idealmodel=ReidIdeal)

    p = 2*101325.0; z = [2.0,5.0];
    q = 0.062744140625
    res_qp3 = qp_flash(model,q,p,z)
    res_qp4 = qp_flash(model,q,p,z./10)
    @test Clapeyron.temperature(res_qp3) ≈ Clapeyron.temperature(res_qp4)

    #VT flash (#331)
    model_a_pr = PR(["water", "oxygen"])
    V_a = 1 # m3
    T = 273.15 + 60 # K
    p_a = 1.2e5 # Pa, 1.2 bar
    n_H2O_a = 1.85e4 # mol H2O
    n_O2_a = 24.08 # mol O2
    sol_fl = vt_flash(model_a_pr, V_a, T, [n_H2O_a, n_O2_a])
    @test V_a ≈ volume(sol_fl)
    water_cpr = cPR(["water"],idealmodel = ReidIdeal)
    @test_throws ArgumentError Clapeyron.VT.speed_of_sound(water_cpr,1e-4,373.15)
    water_cpr_flash = Clapeyron.VT.flash(water_cpr,1e-4,373.15)
    @test_throws ArgumentError speed_of_sound(water_cpr,water_cpr_flash)

    #PH flash with supercritical pure components (#361)
    fluid_model = SingleFluid("Hydrogen")
    T_in = 70               # K
    p_in = 350e5           # Pa
    h_in = enthalpy(fluid_model,p_in,T_in)
    sol_sc = ph_flash(fluid_model,p_in,h_in)
    @test Clapeyron.temperature(sol_sc) ≈ T_in

    #PH Flash where T is in the edge (#373)
    model = cPR(["butane","isopentane"],idealmodel = ReidIdeal)
    p = 101325
    z = [1.0,1.0];
    T = 286.43023797357927 #(0.5*bubble_temperature(model,p,z)[1] + 0.5*dew_temperature(model,p,z)[1])
    h = -50380.604181769755 #Clapeyron.enthalpy(model,p,T,z)
    flash_res_ph = ph_flash(model,p,h,z)
    @test Clapeyron.numphases(flash_res_ph) == 2

    #Inconsistency in flash computations near bubble and dew points (#353)
    fluids =["isopentane","toluene"]
    model = cPR(fluids,idealmodel = ReidIdeal)
    p = 101325
    z = [1.5,1.5]
    T1,T2 = 380, 307.72162335900924 #T1 = 380; T2 = bubble_temperature(model,p,z)[1] - 10
    h1,h2 = 30118.26278687942, -89833.18975112544 #h1 = enthalpy(model,p,T1,z); h2 = enthalpy(model,p,T2,z)
    hrange = range(h1,h2,length=100)
    Trange = similar(hrange)
    for i in eachindex(hrange)
        Ti = Clapeyron.PH.temperature(model,p,hrange[i],z)
        Trange[i] = Ti
        if i > 1
            @test Trange[i] < Trange[i-1] #check that temperature is increasing
            @test isfinite(Ti) #test that there are no NaNs
        end
    end

    #VT flash: water + a tiny amount of hydrogen (#377)
    # content of a cathode separation tank
    n_H2O_c = 0.648e4
    V_c = 0.35
    n_H2_c = 251
    mod_pr = cPR(["water","hydrogen"],idealmodel = ReidIdeal)
    mult_H2 = reverse(0:0.1:5)
    p_tank = similar(mult_H2)
    T_tank = 70 + 273.15
    for (i,mH2) in pairs(mult_H2)
        res_i = vt_flash(mod_pr,V_c,T_tank,[n_H2O_c, exp10(-mH2)*n_H2_c])
        #@test Clapeyron.numphases(res_i) == 2
        #@test pressure(res_i) > 0
        p_tank[i] = pressure(res_i)
    end
    @test count(isnan,p_tank) == 0
    @test issorted(p_tank)

    #394
    fluid394 = cPR(["R134a"],idealmodel=ReidIdeal);
    f394(x) = Clapeyron.PH.temperature(fluid394,101325,x,[1.0]);
    h394 = collect(range(-26617.0,-4282.0,100));
    h394 = -25000.0
    @test iszero(Clapeyron.ForwardDiff.derivative(f394,h394))


    #https://github.com/CoolProp/CoolProp/issues/2622
    model = SingleFluid("R123")
    Mw5 = Clapeyron.molecular_weight(model)
    h5 = 233250.0
    s5 = 1.1049e3
    sm5 = s5*Mw5
    hm5 = h5*Mw5
    p5 = 5e6
    T51 = CoolProp.PropsSI("T","Hmolar",hm5,"P",p5,model)
    T52 = CoolProp.PropsSI("T","H",h5,"P",p5,model)
    T53 = CoolProp.PropsSI("T","Smolar",sm5,"P",p5,model)
    T54 = CoolProp.PropsSI("T","S",s5,"P",p5,model)
    @test T51 == T52
    @test T53 == T54
    @test T53 ≈ 304.88 rtol = 5e-5
    @test T51 ≈ 304.53 rtol = 5e-5

    TUV1 = CoolProp.PropsSI("T","U",29550.0,"D",1000,"water")
    TUV2 = CoolProp.PropsSI("T","U",29550.0,"D",1000,IAPWS95())
    @test TUV1 ≈ TUV2 rtol = 1e-6
    #issue #390
    #=
    model = cPR(["isopentane","toluene"],idealmodel=ReidIdeal)
    z = [0.5,0.5]
    p_crit= 4.1778440598996202e6
    p = collect(range(101325,0.7p_crit,100))
    T_bubble = similar(p)
    T_dew = similar(p)
    s_bubble = similar(p)
    s_dew = similar(p)
    q0 = 0.0
    q1 = 1.0

    for i in eachindex(p)
        res_dew = qp_flash(model,q1,p[i],z)
        T_dew[i] = Clapeyron.temperature(res_dew)
        s_dew[i] = Clapeyron.entropy(model,res_dew)
    end =#
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

    #457
    @test dew_pressure(vdw,T,1)[1] ≈ px rtol = 1e-6
    @test bubble_pressure(vdw,T,1)[1] ≈ px rtol = 1e-6
    @test dew_temperature(vdw,px,1)[1] ≈ T rtol = 1e-6
    @test bubble_temperature(vdw,px,1)[1] ≈ T rtol = 1e-6

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

    #Issue 328
    @test saturation_pressure(cPR("butane"),406.5487245045052)[1] ≈ 2.815259927796967e6 rtol = 1e-6

    #issue 387
    cpr = cPR("Propane",idealmodel = ReidIdeal)
    crit_cpr = crit_pure(cpr)
    @test saturation_temperature(cpr,crit_cpr[2] - 1e3)[1] ≈ 369.88681908031606 rtol = 1e-6
end

@testset "Tproperty/Property" begin
    model1 = cPR(["propane","dodecane"])
    p = 101325.0; T = 300.0;z = [0.5,0.5]
    h_ = enthalpy(model1,p,T,z)
    s_ = entropy(model1,p,T,z)
    @test Tproperty(model1,p,h_,z,enthalpy) ≈ T
    @test Tproperty(model1,p,s_,z,entropy) ≈ T

    model2 = cPR(["propane"])
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

    #issue 409
    fluid409 = cPR(["Propane","R134a"],idealmodel=ReidIdeal);z409 = [1.0,1.0];
    s409 = -104.95768957075641; p409 = 5.910442025416817e6;
    @test Tproperty(fluid409,p409,s409,z409,entropy) ≈ 406.0506318701147 rtol = 1e-6

    model5 = cPR(["R134a","propane"],idealmodel=ReidIdeal)
    @test Clapeyron._Pproperty(model5,450.0,0.03,[0.5,0.5],volume)[2] == :vapour
    @test Clapeyron._Pproperty(model5,450.0,0.03,[0.5,0.5],volume)[2] == :vapour
    @test Clapeyron._Pproperty(model5,450.0,0.00023,[0.5,0.5],volume)[2]  == :eq
    @test Clapeyron._Pproperty(model5,450.0,0.000222,[0.5,0.5],volume)[2]  == :eq
    @test Clapeyron._Pproperty(model5,450.0,0.000222,[0.5,0.5],volume)[2]  == :eq
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
        #not exactly the same results, as activity coefficients are ultimately an approximation of the real helmholtz function.
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

        #413
        fluid413 = cPR(["Propane","Isopentane"],idealmodel=ReidIdeal);
        (p413, y413, method413) = (502277.914581377, [0.9261006181335611, 0.07389938186643885], ChemPotDewTemperature(vol0 = nothing, T0 = nothing, x0 = nothing, noncondensables = nothing, f_limit = 0.0, atol = 1.0e-8, rtol = 1.0e-12, max_iters = 1000, ss = false))
        T413,_,_,_ = Clapeyron.dew_temperature_impl(fluid413,p413,y413,method413)
        @test T413 ≈ 292.1479303719277 rtol = 1e-6
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
        #test if the nonvolatile neq system is being built
        (Tc,vlc,vvc,yc) = bubble_temperature(system2,p,x0,FugBubbleTemperature(itmax_newton = 1, y0 = y0,T0 = T,nonvolatiles = ["decane"]))
        @test Tc isa Number
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
