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


    model = Clapeyron.GERG2008(["nitrogen","methane","ethane","propane","butane","isobutane","pentane"])
    lng_composition = [0.93,92.1,4.64,1.7,0.42,0.32,0.09]
    lng_composition_molar_fractions = lng_composition ./sum(lng_composition)
    @test Clapeyron.molar_density(model,(380.5+101.3)*1000.0,-153.0+273.15,lng_composition_molar_fractions)/1000 ≈ 24.98 rtol = 1E-2
    @test Clapeyron.mass_density(model,(380.5+101.3)*1000.0,-153.0+273.15,lng_composition_molar_fractions) ≈ 440.73 rtol = 1E-2
    @test Clapeyron.molar_density(model,(380.5+101.3)u"kPa",-153.0u"°C",lng_composition_molar_fractions;output=u"mol/L") ≈ 24.98*u"mol/L"  rtol=1E-2
    @test Clapeyron.mass_density(model,(380.5+101.3)u"kPa",-153.0u"°C",lng_composition_molar_fractions;output=u"kg/m^3")  ≈ 440.73*u"kg/m^3" rtol=1E-2
end

@testset "Rachford-Rice" begin
    #error, all Ks > 0, from CalebBell/Chemicals
    zs = [0.0, 0.0885053990596404, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.721469918219507, 0.01742948033685831,
            0.1725952023839942]
    Ks = [192.3625321718047, 105.20070698573475, 76.30532397111489, 37.090890262982, 21.38862102676539, 18.093547012968767,
            10.319129068837443, 9.001962137200403, 4.565198737490148, 340.69153314749224, 269.09234343328467,
            21.858385052861507]
    @test Clapeyron.rr_vle_vapor_fraction(Ks,zs) == Inf

    #here are some tests, from the paper.
    #the paper has errors, z2 has 10 components, 9 and 10 repeated.
    #and after that, they dont include the bmax,bmin in their selected pair.
    z1 = [0.30,0.15, 0.05, 0.02, 0.01, 0.02, 0.02, 0.03, 0.07, 0.33]
    k1 = [3.0E+00, 2.0E+00, 1.1E+00, 8.0E-01, 4.0E-01, 1.0E-01, 5.0E-02, 2.0E-02, 1.0E-02, 1.0E-04]
    idxs1,ks1,zs1 = Clapeyron.rr_find_strongest(k1,z1)
    @test all(in.(idxs1,Ref((1,2,9,10))))
    @test Clapeyron.rr_vle_vapor_fraction(k1,z1) ≈ 0.1769 atol = 1e-4
    z2 = [0.00034825,0.01376300,0.13084000,0.10925000,0.00001000,0.51009000,0.23564000,0.00006000]
    k2 = [5.2663E+02, 5.0400E+01, 1.6463E+00, 8.7450E-01, 1.5589E-01, 3.6588E-02, 2.6625E-02, 4.8918E-06]
    idxs2,ks2,zs2 = Clapeyron.rr_find_strongest(k2,z2)
    #the paper is wrong here, it stablishes (1,2,6,7) ,does not include bmax
    #we still converge to the correct root
    @test all(in.(idxs2,Ref((1,2,6,8))))
    @test Clapeyron.rr_vle_vapor_fraction(k2,z2) ≈ 0.00327 atol = 1e-4

    z3 = [0.0187002, 0.0243002, 0.5419054, 0.0999010, 0.0969010, 0.0400004, 0.0212002, 0.0148001, 0.0741507, 0.0350404, 0.0173602, 0.0157402]
    k3 = [1.32420, 1.12778, 1.22222, 1.11760, 9.88047E-01, 8.94344E-01, 7.87440E-01, 7.43687E-01, 8.11797E-01, 6.93279E-01, 5.09443E-01, 2.28721E-01]
    idxs3,ks3,zs3 = Clapeyron.rr_find_strongest(k3,z3)
    @test all(in.(idxs3,Ref((1,3,9,12))))
    @test Clapeyron.rr_vle_vapor_fraction(k3,z3) ≈ 0.98725 atol = 1e-4

    #testing exact solver,3 comps:
    zs_3 = [0.5, 0.3, 0.2]
    Ks_3 = [1.685, 0.742, 0.532]
    xs_3 = [0.33940869696634357, 0.3650560590371706, 0.2955352439964858]
    ys_3 = [0.5719036543882889, 0.27087159580558057, 0.15722474980613044]
    β_3 = 0.6907302627738544
    β_3_sol = Clapeyron.rr_vle_vapor_fraction(Ks_3,zs_3)
    @test  β_3_sol ≈ β_3
    @test Clapeyron.rr_flash_eval(Ks_3,zs_3,β_3_sol) <= 4*eps(β_3)
    @test Clapeyron.rr_flash_vapor(Ks_3,zs_3,β_3_sol) ≈ ys_3
    @test Clapeyron.rr_flash_liquid(Ks_3,zs_3,β_3_sol) ≈ xs_3

    #Extreme Kvalues, it should give β = 0.99999
    Ks_extreme = [15.464909530837806, 7006.64008090944, 1.8085837711444488, 0.007750676421035811, 30.98450366497431]
    zs_extreme = [0.26562380186293233, 0.04910234829974003, 0.284394553603828, 0.006300023876072552, 0.3945792723574272]
    @test Clapeyron.rr_vle_vapor_fraction(Ks_extreme,zs_extreme) ≈ 0.999999999

    # Case where the evaluated point is right on the boundary
    Ks_eps_0 = [1.2566703532018493e-21, 3.3506275205339295, 1.0300675710905643e-23, 1.706258568414198e-39, 1.6382855298440747e-20]
    zs_eps_0 = [0.13754371891028325, 0.29845155687154623, 0.2546683930289046, 0.08177453852283137, 0.22756179266643456]
    @test abs(Clapeyron.rr_vle_vapor_fraction(Ks_eps_0,zs_eps_0)) < eps(Float64)
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
    #@test Clapeyron.nonzero_extrema(0:3) == (1, 3)
    @test Clapeyron.a_assoc(model_no_comb,V,T,z) ≈ -4.667036481159167 rtol = 1E-6
    @test Clapeyron.a_assoc(model_cr1,V,T,z) ≈ -5.323469194263458 rtol = 1E-6
    @test Clapeyron.a_assoc(model_esd,V,T,z) ≈ -5.323420343872591 rtol = 1E-6
    @test Clapeyron.a_assoc(model_esd_r,V,T,z) ≈ -5.323430326406561 rtol = 1E-6
    @test Clapeyron.a_assoc(model_dufal,V,T,z) ≈ -5.323605338112626 rtol = 1E-6

    #system with strong association:
    fluid = PCSAFT(["water","methanol"]; assoc_options = AssocOptions(combining=:elliott))
    fluid.params.epsilon["water","methanol"] *= (1+0.18)
    v = volume(fluid, 1e5, 160.0, [0.5, 0.5],phase = :l)
    test_assoc_matrix(fluid,v,160.0,[0.5,0.5])
    @test Clapeyron.X(fluid,v,160.0,[0.5,0.5]).v ≈ [0.0011693187791158642, 0.0011693187791158818, 0.0002916842981727242, 0.0002916842981727286] rtol = 1E-8
    #test with bigfloat, we check that all temporary association storage is correctly initialized
    @test Clapeyron.X(fluid,big(v),160.0,[0.5,0.5]).v ≈ [0.0011693187791158642, 0.0011693187791158818, 0.0002916842981727242, 0.0002916842981727286] rtol = 1E-8

    K = [0.0 244.24071691762867 0.0 3.432863727098509; 244.24071691762867 0.0 3.432863727098509 0.0; 0.0 2.2885758180656732 0.0 0.0; 2.2885758180656732 0.0 0.0 0.0]
    @test Clapeyron.assoc_matrix_solve(K) ≈ [0.0562461981664357, 0.0562461981664357, 0.8859564211875895, 0.8859564211875895]
    test_assoc_matrix(K)
    #test some particular 2x2 matrices

    test_assoc_matrix([100.0 200.0; 300.0 400.0])

    test_assoc_matrix([100.0 200.0; 300.0   0.0])
    test_assoc_matrix([100.0   0.0; 200.0 300.0])
    test_assoc_matrix([100.0 200.0; 0.0 300.0])
    test_assoc_matrix([0.0 100.0; 200.0 300.0])

    test_assoc_matrix([0.0 100.0; 200.0 0.0])
    test_assoc_matrix([1000.0 0.0; 0.0 200.0])
    test_assoc_matrix([0.0 1000.0; 0.0 200.0])
    test_assoc_matrix([1000.0 0.0; 200.0 0.0])
    test_assoc_matrix([0.0 0.0; 100.0 200.0])
    test_assoc_matrix([1000.0 3000.0; 0.0 0.0])

    test_assoc_matrix([1000.0 0.0; 0.0 0.0])
    test_assoc_matrix([0.0 1000.0; 0.0 0.0])
    test_assoc_matrix([0.0 0.0; 1000.0 0.0])
    test_assoc_matrix([0.0 0.0; 0.0 1000.0])

    test_assoc_matrix([500.0;;])
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

@testset verbose = true "reference states" begin

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
    @testset "reference states - misc" begin
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
end

@testset verbose = true "Solid Phase Equilibria" begin
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

    #=
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
    end =#
end
GC.gc()

#test for difficult equilibria.
@testset verbose = true "Challenging equilibria" begin

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

        @testset "high conditions - Multifluid" begin
            fluid = MultiFluid(["propane","R134a"])
            #_,pcrit,_ = crit_mix(fluid,[1.0,1.0])
            #p = 0.7*pcrit
            @test dew_temperature(fluid,2.7706485503578815e6,[1.0,1.0])[1] ≈ 335.8362199892292 rtol = 1e-6
        end
    end
    GC.gc()
end

@testset "Partial properties" begin
    model_pem = PR(["hydrogen", "oxygen", "water"])
    z = [0.1,0.1,0.8]
    p,T = 0.95e5,380.15
    for prop in [volume,gibbs_free_energy,helmholtz_free_energy,entropy,enthalpy,internal_energy]
        @test sum(partial_property(model_pem,p,T,z,prop) .* z) ≈ prop(model_pem,p,T,z)
    end
end

@testset "Supercritical lines" begin
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

@testset "Thermodynamic factor" begin 
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
    Γ_ref = thermodynamic_factor_nrtl(act_model, 1e5, T_act, x_act)
    @test size(Γ_act) == (2,2)
    @test all((≈).(Γ_act, Γ_ref, rtol=1e-6)) 
end
