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

    # consistency of the results with/without units
    @test chemical_potential(BasicIdeal(), 1e6u"Pa", 300u"K") == chemical_potential(BasicIdeal(), 1e6, 300)*u"J/mol"
    #@test Clapeyron.x0_psat(model11, 100u"K") == Clapeyron.x0_psat(model11, 100)*u"Pa"
    #@test Clapeyron.x0_sat_pure(model11, 100u"K") == Clapeyron.x0_sat_pure(model11, 100).*(u"m^3",)

    # support for vol0
    modelgergCO2 = GERG2008(["carbon dioxide"])
    @test !isnan(only(Clapeyron.fugacity_coefficient(modelgergCO2, 1u"MPa", 300u"K"; phase=:stable, vol0=0.0023u"m^3")))
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

        alg2 = MichelsenTPFlash(
            equilibrium = :lle, 
            x0 = [0.99999, 0.00001],
            y0 = [0.00001, 0.00009]
        )
        flash2 = tp_flash(system, 101325, 303.15, [0.5, 0.5], alg2)
        act_x2 = activity_coefficient(system, 101325, 303.15, flash2[1][1,:]) .* flash2[1][1,:]
        act_y2 = activity_coefficient(system, 101325, 303.15, flash2[1][2,:]) .* flash2[1][2,:]
        @test Clapeyron.dnorm(act_x2,act_y2) < 1e-8

        #test K0_lle_init initialization
        alg3 = MichelsenTPFlash(
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
        flash4 =  tp_flash(model_vle, 101325, 363.15, [0.5, 0.5], MichelsenTPFlash())
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
    psat_Tcnear = 1.4960621837287119e7 #default solver result
    @test first(Clapeyron.saturation_pressure(model,Tc_near,IsoFugacitySaturation())) ≈ psat_Tcnear rtol = 1e-6
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


    #testset for equilibria bugs
    
end

@testset "Solid Phase Equilibria" begin
    @testset "Pure Solid-Liquid Equilibria" begin
        model = CompositeModel(["methane"]; fluid = SAFTVRMie, solid = SAFTVRSMie)
        
        trp = triple_point(model)
        @test trp[1] ≈ 106.01194395351305 rtol = 1e-6

        sub = sublimation_pressure(model,100.)
        @test sub[1] ≈ 30776.588071307022 rtol = 1e-6

        mel = melting_pressure(model,110.)
        @test mel[1] ≈ 1.126517131058346e7 rtol = 1e-6

        sub = sublimation_temperature(model,1e3)
        @test sub[1] ≈ 78.29626523297529 rtol = 1e-6

        mel = melting_temperature(model,1e5)
        @test mel[1] ≈ 106.02571487518759 rtol = 1e-6

        model2 = CompositeModel("water",solid = SolidHfus, fluid = IAPWS95())
        @test melting_temperature(model2,1e5)[1] ≈ 273.15 rtol = 1e-6
        @test melting_pressure(model2,273.15)[1] ≈ 1e5 rtol = 1e-6
    end

    @testset "Mixture Solid-Liquid Equilibria" begin
        model = CompositeModel([("1-decanol",["CH3"=>1,"CH2"=>9,"OH (P)"=>1]),("thymol",["ACCH3"=>1,"ACH"=>3,"ACOH"=>1,"ACCH"=>1,"CH3"=>2])];liquid=UNIFAC,solid=SolidHfus)
        T = 275.
        p = 1e5
        s1 = sle_solubility(model,p,T,[1.,1.];solute=["1-decanol"])
        s2 = sle_solubility(model,p,T,[1.,1.];solute=["thymol"])
        @test s1[2] ≈ 0.21000625991669147 rtol = 1e-6
        @test s2[2] ≈ 0.3370264930822045 rtol = 1e-6

        (TE,xE) = eutectic_point(model)
        @test TE ≈ 271.97967645045804 rtol = 1e-6
    end
   
    @testset "Solid-Liquid-Liquid Equilibria" begin
        model = CompositeModel(["water","ethanol",("ibuprofen",["ACH"=>4,"ACCH2"=>1,"ACCH"=>1,"CH3"=>3,"COOH"=>1,"CH"=>1])];liquid=UNIFAC,solid=SolidHfus)
        p = 1e5
        T = 323.15
        (s1,s2) = slle_solubility(model,p,T)
        @test s1[3] ≈ 0.0015804179997257882 rtol = 1e-6
    end
end

#test for really really difficult equilibria.
@testset "challenging equilibria" begin
       
    #see https://github.com/ClapeyronThermo/Clapeyron.jl/issues/173
    @testset "VTPR - 1" begin
        #=
        carbon monoxide is supercritical.
        =#

        system = VTPR(["carbon monoxide","carbon dioxide"])
        @test_broken Clapeyron.bubble_pressure(system,218.15,[1e-5,1-1e-5])[1] ≈ 1.1373024916997014e6 rtol = 1e-4
    end



    #see https://github.com/ClapeyronThermo/Clapeyron.jl/issues/172
    @testset "PCSAFT - 1" begin
        #=
        really near critical temperature of the mixture
        seems that was fixed by passing the initial point to the x0_bubble_pressure function
        =#
        x = [0.96611,0.01475,0.01527,0.00385]
        T = 202.694
        v0 = [-4.136285855713797, -4.131888756537859, 0.9673991775701574, 0.014192499147585259, 0.014746430039492817, 0.003661893242764558]
        model = PCSAFT(["methane","butane","isobutane","pentane"])
        # @test_broken bubble_pressure(model,T,x;v0 = v0)[1] ≈ 5.913118531569793e6 rtol = 1e-4
        # FIXME: The test does not yield the same value depending on the OS and the julia version
    end

    @testset "saturation pressures" begin
        model1 = PCSAFT("water")
        Tc1,_,_ = crit_pure(model1)
        T1 = 0.995Tc1
        @test Clapeyron.saturation_pressure(model1,T1,crit_retry = false)[1] ≈ 3.542008160105954e7 rtol = 1e-6
        
        model2 = PCSAFT("eicosane")
        Tc2,_,_ = crit_pure(model2)
        T2 = 0.995Tc2
        if Base.VERSION >= v"1.7" #this test fails on mac, julia 1.6
            @test Clapeyron.saturation_pressure(model2,T2,crit_retry = false)[1] ≈ 1.3931662325210017e6 rtol = 1e-6
        end

        #https://github.com/ClapeyronThermo/Clapeyron.jl/issues/237
        model3 = SAFTVRMie("heptacosane",userlocations = (Mw = 380.44,segment = 2.0,sigma = 3.0,lambda_a = 6.0,lambda_r = 20.01,epsilon = 200.51))
        @test Clapeyron.saturation_pressure(model3,94.33,crit_retry = false)[1] ≈ 2.8668634416924506 rtol = 1e-6
    end

end

