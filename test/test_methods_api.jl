@testset "Unitful Methods" begin
    model11 = GERG2008(["methane"])
    model10 = GERG2008(["butane"])
    #example 3.11 abott van ness, 7th ed.
    #pressure. 189 atm with CS compressibility relation
    p11 = 185.95465583962599u"atm"
    v11 = 2u"ft^3"
    T11 = 122u"°F"
    n11 = 453.59237u"mol" #1 lb-mol
    z11 = 0.8755772456569365 #t0.89 from CS compressibility relation
    @test Clapeyron.pressure(model11,v11,T11,n11,output = u"atm") ≈ p11
    @test Clapeyron.pressure(model11,v11,T11,[n11],output = u"atm") ≈ p11
    @test Clapeyron.compressibility_factor(model11,v11,T11,n11) ≈ z11 rtol = 1E-6
    @test Clapeyron.compressibility_factor(model11,p11,T11,n11) ≈ z11 rtol = 1E-6

    #example 3.10 abott van ness, 7th ed.
    #volume, 1480 cm3, with CS virial correlation
    p10 = 25u"bar"
    T10 = 510u"K"
    Tc10 = 425.75874890467253u"K"
    pc10 = 3.830319495176967e6u"Pa"
    R = (Clapeyron.R̄)u"J/(K*mol)"
    v10 = 1478.2681326033257u"cm^3"
    @test volume(model10,p10,T10,output=u"cm^3") ≈ v10 rtol = 1E-6
    #generalized pitzer CS virial gives -0.220
    @test Clapeyron.second_virial_coefficient(model10,T10)*pc10/(R*Tc10) |> Unitful.ustrip ≈ -0.22346581496303466 rtol = 1E-6

    #example 3.13, abbott and van ness, 7th ed.
    model13 = PR(["ammonia"],translation = RackettTranslation)
    v13 = 26.545208120801895u"cm^3"
    T13 = 310u"K"
    #experimental value is 29.14 cm3/mol. PR default is ≈ 32, Racckett overcorrects
    @test saturation_pressure(model13,T13,output = (u"atm",u"cm^3",u"cm^3"))[2] ≈ v13 rtol = 1E-6
    @test Clapeyron.pip(model13,v13,T13) > 1 #check if is a liquid phase

    #problem 3.1 abbott and van ness, 7th ed.
    model31 = IAPWS95()
    v31 = volume(model31,1u"bar",50u"°C")
    #experimental value is 44.18e-6. close enough.

    @test isothermal_compressibility(model31,1u"bar",50u"°C",output = u"bar^-1") ≈ 44.17306906730427e-6u"bar^-1" rtol = 1E-6
    @test isothermal_compressibility(model31,1u"bar",50u"°C",output = u"bar^-1") ≈ 44.17306906730427e-6u"bar^-1" rtol = 1E-6
    #enthalpy of vaporization of water at 100 °C
    @test enthalpy_vap(model31,100u"°C",output = u"kJ") ≈ 40.64971775824767u"kJ" rtol = 1E-6
end

@testset "association" begin
    no_comb_sparse = Clapeyron.AssocOptions(combining = :nocombining, dense = false)
    no_comb_dense = Clapeyron.AssocOptions(combining = :nocombining, dense = true)
    esd = Clapeyron.AssocOptions(combining = :esd)
    esd_r = Clapeyron.AssocOptions(combining = :elliott_runtime)
    cr1 = Clapeyron.AssocOptions(combining = :cr1)
    
    model_no_comb_sparse = PCSAFT(["methanol","ethanol"],assoc_options = no_comb_sparse)
    model_no_comb_dense = PCSAFT(["methanol","ethanol"],assoc_options = no_comb_dense)
    model_cr1 = PCSAFT(["methanol","ethanol"],assoc_options = cr1)
    model_esd = PCSAFT(["methanol","ethanol"],assoc_options = esd)
    model_esd_r = PCSAFT(["methanol","ethanol"],assoc_options = esd_r)

    V = 5e-5
    T = 298.15
    z = [0.5,0.5]
    @test Clapeyron.nonzero_extrema(0:3) == (1, 3)
    @test Clapeyron.a_assoc(model_no_comb_sparse,V,T,z) ≈ -4.667036481159167  rtol = 1E-6
    @test Clapeyron.a_assoc(model_no_comb_sparse,V,T,z) ≈ Clapeyron.a_assoc(model_no_comb_dense,V,T,z)  rtol = 1E-6
    @test Clapeyron.a_assoc(model_cr1,V,T,z) ≈ -5.323469194263458  rtol = 1E-6
    @test Clapeyron.a_assoc(model_esd,V,T,z) ≈ -5.323420343872591  rtol = 1E-6
    @test Clapeyron.a_assoc(model_esd_r,V,T,z) ≈ -5.323430326406561  rtol = 1E-6
end

@testset "tpd" begin
    system = PCSAFT(["water","cyclohexane"])
    T = 298.15
    p = 1e5
    phases,tpds,symz,symw = Clapeyron.tpd(system,p,T,[0.5,0.5])
    @test tpds[1] ≈ -0.8370113547074933  rtol = 1e-6
    GC.gc()
end


@testset "Tp flash algorithms" begin
    #this is a VLLE equilibria
    system = PCSAFT(["water","cyclohexane","propane"])
    T = 298.15
    p = 1e5
    z = [0.333, 0.333, 0.334]

    @testset "RR Algorithm" begin
        method = RRTPFlash()
        @test Clapeyron.tp_flash(system, p, T, z, method)[3] ≈ -6.539976318817461 rtol = 1e-6
    end
    GC.gc()

    @testset "DE Algorithm" begin
        method = DETPFlash(numphases=3)
        @test Clapeyron.tp_flash(system, p, T, z, method)[3] ≈ -6.759674475174073 rtol = 1e-6
    end
    GC.gc()

    @testset "Michelsen Algorithm" begin

        x0 = [0.9997755902156433, 0.0002244097843566859, 0.0]
        y0 = [6.425238373915699e-6, 0.9999935747616262, 0.0]
        method = MichelsenTPFlash(x0 = x0, y0 = y0, equilibrium= :lle)
        @test Clapeyron.tp_flash(system, p, T, [0.5,0.5,0.0],method)[3] ≈ -7.577270350886795 rtol = 1e-6

        method2 = MichelsenTPFlash(x0 = x0, y0 = y0, equilibrium = :lle, ss_iters = 1, second_order = false)
        @test Clapeyron.tp_flash(system, p, T, [0.5,0.5,0.0],method2)[3] ≈ -7.577270350886795 rtol = 1e-6
    
        method3 = MichelsenTPFlash(x0 = x0, y0 = y0, equilibrium = :lle, ss_iters = 1,second_order = true)
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
end


@testset "Saturation Methods" begin
    model = PR(["water"])
    vdw = vdW(["water"])
    p0 = 1e5
    T = 373.15
    p,vl,vv = Clapeyron.saturation_pressure(model,T) #default
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
    psat_Tcnear = 1.4960621837287119e7 #default solver result
    @test first(Clapeyron.saturation_pressure(model,Tc_near,IsoFugacitySaturation())) ≈ psat_Tcnear rtol = 1e-6
    #Test that IsoFugacity fails over critical point
    @test isnan(first(Clapeyron.saturation_pressure(model,1.1*647.096,IsoFugacitySaturation())))
    GC.gc()

    #SuperAncSaturation
    p5,vl5,vv5 = Clapeyron.saturation_pressure_impl(model,T,SuperAncSaturation())
    @test p5 ≈ p rtol = 1e-6
    @test @inferred Clapeyron.saturation_pressure_impl(vdw,T,SuperAncSaturation())[1] ≈ px
    GC.gc()

    #AntoineSat
    @test Clapeyron.saturation_temperature(model,p0,AntoineSaturation(T0 = 400.0))[1] ≈ 374.2401401001685 rtol = 1e-6
    @test Clapeyron.saturation_temperature(model,p0,AntoineSaturation(vl = vl5,vv = vv5))[1] ≈ 374.2401401001685 rtol = 1e-6
    @test_throws Any Clapeyron.saturation_temperature(model,p0,AntoineSaturation(vl = vl5,T0 = 400))
    GC.gc()

    #ClapeyronSat
    @test Clapeyron.saturation_temperature(model,p0,ClapeyronSaturation())[1] ≈ 374.2401401001685 rtol = 1e-6
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
        @test Clapeyron.dew_temperature(system1,p2,z,Clapeyron.FugDewTemperature(itmax_newton = 1))[1] ≈ Tres2 rtol = 1E-6
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
