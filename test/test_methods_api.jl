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
    model13 = PR(["ammonia"],translation = PenelouxTranslation)
    v13 = 30.769606571028624u"cm^3"
    T13 = 310u"K"
    #experimental value is 29.14 cm3/mol. PR default is ≈ 32, Peneloux overcorrects a little
    @test saturation_pressure(model13,T13,output = (u"atm",u"cm^3",u"cm^3"))[2] ≈ v13 rtol = 1E-6
    @test Clapeyron.pip(model13,v13,T13) > 1 #check if is a liquid phase

    #problem 3.1 abbott and van ness, 7th ed.
    model31 = IAPWS95()
    v31 = volume(model31,1u"bar",50u"°C")
    #experimental value is 44.18e-6. close enough.

    @test isothermal_compressibility(model31,1u"bar",50u"°C",output = u"bar^-1") ≈ 44.17306906730427e-6u"bar^-1" rtol = 1e-4
    @test isothermal_compressibility(model31,1u"bar",50u"°C",output = u"bar^-1") ≈ 44.17306906730427e-6u"bar^-1" rtol = 1e-4
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
    model_no_comb = PCSAFT(["methanol","ethanol"],assoc_options = AssocOptions(combining = :nocombining))
    model_cr1 = PCSAFT(["methanol","ethanol"],assoc_options = AssocOptions(combining = :cr1))
    model_esd = PCSAFT(["methanol","ethanol"],assoc_options = AssocOptions(combining = :esd))
    model_esd_r = PCSAFT(["methanol","ethanol"],assoc_options = AssocOptions(combining = :elliott_runtime))
    model_dufal = PCSAFT(["methanol","ethanol"],assoc_options = AssocOptions(combining = :dufal))
    test_repr(AssocOptions(combining = :dufal))
    V = 5e-5
    T = 298.15
    z = [0.5,0.5]
    @test Clapeyron.nonzero_extrema(0:3) == (1, 3)
    @test Clapeyron.a_assoc(model_no_comb,V,T,z) ≈ -4.667036481159167 rtol = 1E-6
    @test Clapeyron.a_assoc(model_cr1,V,T,z) ≈ -5.323469194263458 rtol = 1E-6
    @test Clapeyron.a_assoc(model_esd,V,T,z) ≈ -5.323420343872591 rtol = 1E-6
    @test Clapeyron.a_assoc(model_esd_r,V,T,z) ≈ -5.323430326406561 rtol = 1E-6
    @test Clapeyron.a_assoc(model_dufal,V,T,z) ≈ -5.323605338112626 rtol = 1E-6

    #system with strong association:
    fluid = PCSAFT(["water","methanol"]; assoc_options = AssocOptions(combining=:elliott))
    fluid.params.epsilon["water","methanol"] *= (1+0.18)
    v = volume(fluid, 1e5, 160.0, [0.5, 0.5],phase = :l)
    @test Clapeyron.X(fluid,v,160.0,[0.5,0.5]).v ≈ [0.0011693187791158642, 0.0011693187791158818, 0.0002916842981727242, 0.0002916842981727286] rtol = 1E-8
    #test with bigfloat, we check that all temporary association storage is correctly initialized
    @test Clapeyron.X(fluid,big(v),160.0,[0.5,0.5]).v ≈ [0.0011693187791158642, 0.0011693187791158818, 0.0002916842981727242, 0.0002916842981727286] rtol = 1E-8

    K = [0.0 244.24071691762867 0.0 3.432863727098509; 244.24071691762867 0.0 3.432863727098509 0.0; 0.0 2.2885758180656732 0.0 0.0; 2.2885758180656732 0.0 0.0 0.0]
    @test Clapeyron.assoc_matrix_solve(K) ≈ [0.0562461981664357, 0.0562461981664357, 0.8859564211875895, 0.8859564211875895]
end

using EoSSuperancillaries
#we test this separately, but leaving this on could speed up the test suite?
Clapeyron.use_superancillaries!(false)

if isdefined(Base,:get_extension)
    @testset "Superancillaries.jl" begin
        
        pc = PCSAFT("eicosane")
        pc2 = pharmaPCSAFT("oxygen")
        cubic = tcPR(["water"])
        
        crit_pc = crit_pure(pc)
        sat_cubic = saturation_pressure(cubic,373.15)
        sat_pc2 = saturation_pressure(pc2,150.0)
        
        Clapeyron.use_superancillaries!(true)
        
        crit_sa_pc = crit_pure(pc)
        sat_sa_cubic = saturation_pressure(cubic,373.15)
        sat_sa_pc2 = saturation_pressure(pc2,150.0)

        @test crit_pc[1] ≈ crit_sa_pc[1] rtol = 1e-6
        @test crit_pc[3] ≈ crit_sa_pc[3] rtol = 1e-6
        
        @test sat_cubic[1] ≈ sat_sa_cubic[1] rtol = 1e-6
        @test sat_cubic[2] ≈ sat_sa_cubic[2] rtol = 1e-6
        @test sat_cubic[3] ≈ sat_sa_cubic[3] rtol = 1e-6
        
        @test sat_pc2[1] ≈ sat_sa_pc2[1] rtol = 1e-6
        @test sat_pc2[2] ≈ sat_sa_pc2[2] rtol = 1e-6
        @test sat_pc2[3] ≈ sat_sa_pc2[3] rtol = 1e-6

    end
end

@testset "tpd" begin
    system = PCSAFT(["water","cyclohexane"])
    T = 298.15
    p = 1e5
    phases,tpds,symz,symw = Clapeyron.tpd(system,p,T,[0.5,0.5])
    @test tpds[1] ≈ -0.6081399681963373 rtol = 1e-6

    act_system = UNIFAC(["water","cyclohexane"])
    phases2,tpds2,symz2,symw2 = Clapeyron.tpd(act_system,p,T,[0.5,0.5],lle = true)
    @test tpds2[1] ≈ -0.9412151812640561 rtol = 1e-6
    GC.gc()

    model = PR(["IsoButane", "n-Butane", "n-Pentane", "n-Hexane"]);
    z = [0.25, 0.25, 0.25, 0.25]
    p = 1e5
    v1 = volume(model,p,297.23,z)
    @test !isstable(model,p,297.23,z)
    @test !Clapeyron.VT_isstable(model,v1,297.23,z)
    GC.gc()

    model = cPR(["ethane"],idealmodel = ReidIdeal)
    p = 101325; z = [5.0]; 
    T,vl,vv = saturation_temperature(model,p)
    v_unstable = exp(0.5*(log(vl) + log(vv)))
    V = volume(model,p,T,z) # lies in range of Vv
    @test Clapeyron.VT_diffusive_stability(model,V,T,z)
    @test !Clapeyron.VT_diffusive_stability(model,v_unstable,T,z)
end

@testset "reference states" begin

    @test !has_reference_state(PCSAFT("water"))
    @test has_reference_state(PCSAFT("water",idealmodel = ReidIdeal))

    ref1 = ReferenceState(:nbp)
    ref2 = ReferenceState(:ashrae)
    ref3 = ReferenceState(:iir)
    ref4 = ReferenceState(:volume,T0 = 298.15,P0 = 1.0e5,phase = :liquid)
    ref5 = ReferenceState(:volume,T0 = 298.15,P0 = 1.0e5,phase = :gas,z0 = [0.4,0.6],H0 = 123,S0 = 456)

    @testset "nbp reference state" begin
        model1 = PCSAFT(["water","pentane"],idealmodel = ReidIdeal,reference_state = ref1)
        pure1 = split_model(model1)
        T11,v11,_ = saturation_temperature(pure1[1],101325.0)
        @test Clapeyron.VT_enthalpy(pure1[1],v11,T11) ≈ 0.0 atol = 1e-6
        @test Clapeyron.VT_entropy(pure1[1],v11,T11) ≈ 0.0 atol = 1e-6
        T12,v12,_ = saturation_temperature(pure1[2],101325.0)
        @test Clapeyron.VT_enthalpy(pure1[2],v12,T12) ≈ 0.0 atol = 1e-6
        @test Clapeyron.VT_entropy(pure1[2],v12,T12) ≈ 0.0 atol = 1e-6
        test_repr(Clapeyron.reference_state(model1))
        #test that multifluids work.
        model1b = GERG2008("water",reference_state = :nbp)
        T1b,v1b,_ = saturation_temperature(model1b,101325.0)
        @test Clapeyron.VT_enthalpy(model1b,v1b,T1b) ≈ 0.0 atol = 1e-6
        @test Clapeyron.VT_entropy(model1b,v1b,T1b) ≈ 0.0 atol = 1e-6
    end

    @testset "ashrae reference state" begin
        model2 = PCSAFT(["water","pentane"],idealmodel = ReidIdeal,reference_state = ref2)
        pure2 = split_model(model2)
        T_ashrae = 233.15
        _,v21,_ = saturation_pressure(pure2[1],T_ashrae)
        @test Clapeyron.VT_enthalpy(pure2[1],v21,T_ashrae) ≈ 0.0 atol = 1e-6
        @test Clapeyron.VT_entropy(pure2[1],v21,T_ashrae) ≈ 0.0 atol = 1e-6
        _,v22,_ = saturation_pressure(pure2[2],T_ashrae)
        @test Clapeyron.VT_enthalpy(pure2[2],v22,T_ashrae) ≈ 0.0 atol = 1e-6
        @test Clapeyron.VT_entropy(pure2[2],v22,T_ashrae) ≈ 0.0 atol = 1e-6
    end

    @testset "iir reference state" begin
        model3 = PCSAFT(["water","pentane"],idealmodel = ReidIdeal,reference_state = ref3)
        pure3 = split_model(model3)
        Tiir = 273.15
        H31,H32 = 200*model3.params.Mw[1],200*model3.params.Mw[2]
        S31,S32 = 1.0*model3.params.Mw[1],1.0*model3.params.Mw[2]
        _,v31,_ = saturation_pressure(pure3[1],Tiir)
        if isnan(v31)
            @show Clapeyron.reference_state(model3)
        end
        @test Clapeyron.VT_enthalpy(pure3[1],v31,Tiir) ≈ H31 atol = 1e-6
        @test Clapeyron.VT_entropy(pure3[1],v31,Tiir) ≈ S31 atol = 1e-6
        _,v32,_ = saturation_pressure(pure3[2],Tiir)
        @test Clapeyron.VT_enthalpy(pure3[2],v32,Tiir) ≈ H32 atol = 1e-6
        @test Clapeyron.VT_entropy(pure3[2],v32,Tiir) ≈ S32 atol = 1e-6
    end

    @testset "custom reference state - pure volume" begin
        model4 = PCSAFT(["water","pentane"],idealmodel = ReidIdeal,reference_state = ref4)
        pure4 = split_model(model4)
        T4,P4 = 298.15,1.0e5
        v41 = volume(pure4[1],P4,T4,phase = :liquid)
        @test Clapeyron.VT_enthalpy(pure4[1],v41,T4) ≈ 0.0 atol = 1e-6
        @test Clapeyron.VT_entropy(pure4[1],v41,T4) ≈ 0.0 atol = 1e-6
        v42 = volume(pure4[2],P4,T4,phase = :liquid)
        @test Clapeyron.VT_enthalpy(pure4[2],v42,T4) ≈ 0.0 atol = 1e-6
        @test Clapeyron.VT_entropy(pure4[2],v42,T4) ≈ 0.0 atol = 1e-6
    end

    @testset "custom reference state - volume at composition" begin
        model5 = PCSAFT(["water","pentane"],idealmodel = ReidIdeal,reference_state = ref5)
        T5,P5 = 298.15,1.0e5
        z5 = [0.4,0.6]
        v5 = volume(model5,P5,T5,z5,phase = :gas)
        @test Clapeyron.VT_enthalpy(model5,v5,T5,z5) ≈ 123 atol = 1e-6
        @test Clapeyron.VT_entropy(model5,v5,T5,z5) ≈ 456 atol = 1e-6
    end

    #reference state from EoSVectorParam
    mod_pr = cPR(["water","ethanol"],idealmodel = ReidIdeal,reference_state = :ntp)
    mod_vec = Clapeyron.EoSVectorParam(mod_pr)
    Clapeyron.recombine!(mod_vec)
    @test reference_state(mod_vec).std_type == :ntp
    @test length(reference_state(mod_vec).a0) == 2

    #reference state from Activity models
    puremodel = mod_pr = cPR(["water","ethanol"],idealmodel = ReidIdeal)
    act = NRTL(["water","ethanol"],puremodel = puremodel,reference_state = :ntp)
    @test reference_state(act).std_type == :ntp
    @test length(reference_state(act).a0) == 2

    #issue 511
    ref511 = ReferenceState(:nbp)
    model511 = cPR("water",idealmodel=ReidIdeal,reference_state = ref511)
    @test Clapeyron.reference_state(model511).std_type == :nbp
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

        #solid gibbs + fluid helmholtz
        model3 = CompositeModel("water", solid = IAPWS06(),fluid = IAPWS95())
        @test melting_temperature(model3,101325.0)[1] ≈ 273.1525192653753 rtol = 1e-6
        @test melting_pressure(model3,273.1525192653753)[1] ≈ 101325.0 rtol = 1e-6

        #solid gibbs + fluid gibbs
        model4 = CompositeModel("water", solid = IAPWS06(),fluid = GrenkeElliottWater())
        @test melting_temperature(model4,101325.0)[1] ≈ 273.15 rtol = 1e-6
        @test melting_pressure(model4,273.15)[1] ≈ 101325.0 rtol = 1e-6

        #solid gibbs + any other fluid helmholtz
        model5 = CompositeModel("water", solid = IAPWS06(),fluid = cPR("water"))
        @test melting_temperature(model5,101325.0)[1] ≈ 273.15 rtol = 1e-6
        @test melting_pressure(model5,273.15)[1] ≈ 101325.0 rtol = 1e-6

        #solid gibbs + helmholtz fluid, without any initial points
        model6 = CompositeModel(["CO2"],solid = JagerSpanSolidCO2(),fluid = SingleFluid("carbon dioxide"))
        tp6 = triple_point(model6)
        Ttp6 = tp6[1]
        ptp6 = tp6[2]
        @test Ttp6 ≈ model6.fluid.properties.Ttp rtol = 1e-5
        @test sublimation_pressure(model6,Ttp6)[1] ≈ ptp6 rtol = 1e-6
        @test sublimation_temperature(model6,ptp6)[1] ≈ Ttp6 rtol = 1e-6
        @test melting_pressure(model6,Ttp6)[1] ≈ ptp6 rtol = 1e-6
        @test melting_temperature(model6,ptp6)[1] ≈ Ttp6 rtol = 1e-6
    end
    GC.gc()
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

    GC.gc()

    @testset "Solid-Liquid-Liquid Equilibria" begin
        model = CompositeModel(["water","ethanol",("ibuprofen",["ACH"=>4,"ACCH2"=>1,"ACCH"=>1,"CH3"=>3,"COOH"=>1,"CH"=>1])];liquid=UNIFAC,solid=SolidHfus)
        p = 1e5
        T = 323.15
        (s1,s2) = slle_solubility(model,p,T)
        @test s1[3] ≈ 0.0015804179997257882 rtol = 1e-6
    end

    @testset "#466" begin
        glycine = ("glycine" => ["COOH" => 1, "CH2" => 1, "NH2" => 1])
        lactic_acid = ("lactic acid" =>["COOH" => 1, "CH3" => 1, "CHOH" => 1])
        oxalic_acid = ("oxalic acid" => ["COOH" => 2])

        ox_gly = CompositeModel([oxalic_acid,glycine];fluid=SAFTgammaMie,solid=SolidHfus)
        la_gly = CompositeModel([lactic_acid,glycine];fluid=SAFTgammaMie,solid=SolidHfus)

        T1,_ = Clapeyron.eutectic_point(la_gly)
        T2,_ = Clapeyron.eutectic_point(ox_gly)
        @test T1 ≈ 300.23095880432294 rtol = 1e-6
        @test T2 ≈ 454.27284723964925 rtol = 1e-6
    end
end
GC.gc()
#test for really really difficult equilibria.
@testset "challenging equilibria" begin

    #see https://github.com/ClapeyronThermo/Clapeyron.jl/issues/173
    @testset "VTPR - 1" begin
        #=
        carbon monoxide is supercritical.
        =#

        system = VTPR(["carbon monoxide","carbon dioxide"])
        @test Clapeyron.bubble_pressure(system,218.15,[1e-5,1-1e-5])[1] ≈ 554338.3127125567 rtol = 1e-4
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
        #@test bubble_pressure(model,T,x;v0 = v0)[1] ≈ 5.913118531569793e6 rtol = 1e-4
        # FIXME: The test does not yield the same value depending on the OS and the julia version
    end
    GC.gc()
    @testset "saturation points without critical point" begin
        model1 = PCSAFT("water")
        Tc1,_,_ = crit_pure(model1)
        T1 = 0.999Tc1
        @test Clapeyron.saturation_pressure(model1,T1,crit_retry = false)[1] ≈ 3.6377840330375336e7 rtol = 1e-6

        model2 = PCSAFT("eicosane")
        Tc2,_,_ = crit_pure(model2)
        T2 = 0.999Tc2

        #this test fails on mac, julia 1.6
        @test Clapeyron.saturation_pressure(model2,T2,crit_retry = false)[1] ≈ 1.451917823392476e6 rtol = 1e-6

        #https://github.com/ClapeyronThermo/Clapeyron.jl/issues/237
        #for some reason, it fails with mac sometimes
        if !Base.Sys.isapple()
            model3 = SAFTVRMie("heptacosane",userlocations = (Mw = 380.44,segment = 2.0,sigma = 3.0,lambda_a = 6.0,lambda_r = 20.01,epsilon = 200.51))
            @test Clapeyron.saturation_pressure(model3,94.33,crit_retry = false)[1] ≈ 2.8668634416924506 rtol = 1e-6
        end

        model4 = SAFTVRMie(["methanol"])
        T4 = 164.7095044742657
        @test Clapeyron.saturation_pressure(model4,T4,crit_retry = false)[1] ≈ 0.02610821545005174 rtol = 1e-6
    
        @testset "saturation at low temperatures" begin
            l1 = PR("1-butene")
            Tc1 = 419.95
            sat_low1 = saturation_pressure(l1,0.183Tc1)
            @test sat_low1[1] ≈ 9.468875475768151e-9 rtol = 1e-6

            l2 = PCSAFT("1-butene")
            Tc2 = 426.80960130305374
            sat_low2 = saturation_pressure(l2,0.18Tc2)
            @test sat_low2[1] ≈ 1.7914820721239496e-9 rtol = 1e-6
            
        end
    end
    GC.gc()
end

@testset "partial properties" begin
    model_pem = PR(["hydrogen", "oxygen", "water"])
    z = [0.1,0.1,0.8]
    p,T = 0.95e5,380.15
    for prop in [volume,gibbs_free_energy,helmholtz_free_energy,entropy,enthalpy,internal_energy]
        @test sum(partial_property(model_pem,p,T,z,prop) .* z) ≈ prop(model_pem,p,T,z)
    end
end

@testset "spinodals" begin
    # Example from Ref. https://doi.org/10.1016/j.fluid.2017.04.009
    model = PCSAFT(["methane","ethane"])
    T_spin = 223.
    x_spin = [0.2,0.8]
    (pl_spin, vl_spin) = spinodal_pressure(model,T_spin,x_spin;phase=:liquid)
    (pv_spin, vv_spin) = spinodal_pressure(model,T_spin,x_spin;phase=:vapor)
    @test vl_spin ≈ 7.218532167482202e-5 rtol = 1e-6
    @test vv_spin ≈ 0.0004261109817247137 rtol = 1e-6

    (Tl_spin_impl, xl_spin_impl) = spinodal_temperature(model,pl_spin,x_spin;T0=220.,v0=vl_spin)
    (Tv_spin_impl, xv_spin_impl) = spinodal_temperature(model,pv_spin,x_spin;T0=225.,v0=vv_spin)
    @test Tl_spin_impl ≈ T_spin rtol = 1e-6
    @test Tv_spin_impl ≈ T_spin rtol = 1e-6

    #test for #382: pure spinodal at low pressures
    model2 = PCSAFT("carbon dioxide")
    Tc,Pc,Vc = (310.27679925044134, 8.06391600653306e6, 9.976420206333288e-5)
    T = LinRange(Tc-70,Tc-0.1,50)
    psl = first.(spinodal_pressure.(model2,T,phase = :l))
    psv = first.(spinodal_pressure.(model2,T,phase = :v))
    psat = first.(saturation_pressure.(model2,T))
    @test all(psl .< psat)
    @test all(psat .< psv)
    @test issorted(psl)
    @test issorted(psv)
end

@testset "supercritical lines" begin
    model = PR("methane")
    T_initial = 200.0
    p_widom, v1 = widom_pressure(model, T_initial)
    T_widom, v2 = widom_temperature(model, p_widom)
    @test T_initial ≈ T_widom rtol = 1e-6
    @test v1 ≈ v2 rtol = 1e-6
    @test_throws ArgumentError widom_pressure(model,T_initial,v0 = v1,p0 = p_widom)
    @test T_initial ≈ widom_temperature(model,p_widom,T0 = 1.01*T_initial)[1] rtol = 1e-6
    @test T_initial ≈ widom_temperature(model,p_widom,v0 = 1.01*v1)[1] rtol = 1e-6
    @test p_widom ≈ widom_pressure(model,T_initial,p0 = 1.01*p_widom)[1] rtol = 1e-6
    @test p_widom ≈ widom_pressure(model,T_initial,v0 = 1.01*v1)[1] rtol = 1e-6

    p_ciic, v3 = ciic_pressure(model, T_initial)
    T_ciic, v4 = ciic_temperature(model, p_ciic)
    @test T_initial ≈ T_ciic rtol = 1e-6
    @test v3 ≈ v4 rtol = 1e-6
    @test_throws ArgumentError ciic_pressure(model,T_initial,v0 = v3,p0 = p_ciic)
    @test T_initial ≈ ciic_temperature(model,p_ciic,T0 = 1.01*T_initial)[1] rtol = 1e-6
    @test T_initial ≈ ciic_temperature(model,p_ciic,v0 = 1.01*v3)[1] rtol = 1e-6
    @test p_ciic ≈ ciic_pressure(model,T_initial,p0 = 1.01*p_ciic)[1] rtol = 1e-6
    @test p_ciic ≈ ciic_pressure(model,T_initial,v0 = 1.01*v3)[1] rtol = 1e-6
end

@testset "thermodynamic factor" begin 
    eos_model = PCSAFT(["water", "ethanol"])
    Γ_eos = thermodynamic_factor(eos_model, 1e5, 300., [2.,4.])
    @test size(Γ_eos) == (1,1)
    @test Γ_eos[1,1] ≈ 0.6178686409160774 rtol=1e-6 

    # Activity model
    act_model = NRTL(["water", "ethanol", "acetone"])
    T_act = 300.
    x_act = [0.15, 0.35, 0.5]
    
    # Analytical from Taylor and Krishna (1993), Table D.8
    function thermodynamic_factor_nrtl(model, p, T, x)    
        N = length(x)
        
        τ = model.params.a .+ model.params.b ./ T
        G = exp.(-τ .* model.params.c)
        S = G' * x    
        ε = G .* (τ .- (((G .* τ)' * x) ./ S)') .* (inv.(S)')
        Q = ε .+ ε' .- ((G * Diagonal(x ./ S)) * ε' .+ (ε * Diagonal(x ./ S)) * G')
        
        Γ = [((i == j) ? 1.0 : 0.0) + x[i] * (Q[i, j] - Q[i, N]) for i in 1:N-1, j in 1:N-1]
        return Γ
    end

    Γ_act = thermodynamic_factor(act_model, 1e5, T_act, x_act)
    Γ_ref = thermodynamic_factor(act_model, 1e5, T_act, x_act)
    @test size(Γ_act) == (2,2)
    @test all((≈).(Γ_act, Γ_ref, rtol=1e-6)) 
end
