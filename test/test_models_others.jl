

@testset "Activity models" begin
    @printline
    let T = 333.15, V = 1e-3, p = 1e5, z = [0.5, 0.5], z1 = Clapeyron.SA[1.0], z2 = [0.5,0.5], z3 = [0.333, 0.333, 0.333];
    @testset "Margules" begin
        system = Margules(["methanol","water"])
       #= @test Clapeyron.activity_coefficient(system,p,T,z)[1] #≈ 1.530046633499114 rtol = 1e-6
        @test Clapeyron.activity_coefficient(system,p,T,z) #≈ Clapeyron.test_activity_coefficient(system,p,T,z)  rtol = 1e-6
        =#@test Clapeyron.excess_gibbs_free_energy(system,p,T,z) ≈ Clapeyron.test_excess_gibbs_free_energy(system,p,T,z)  rtol = 1e-6
    end
    @testset "VanLaar" begin
        system = VanLaar(["methanol","water"])
     #=   @test Clapeyron.activity_coefficient(system,p,T,z)[1] #≈ 1.530046633499114 rtol = 1e-6
        @test Clapeyron.activity_coefficient(system,p,T,z) #≈ Clapeyron.test_activity_coefficient(system,p,T,z)  rtol = 1e-6
        =#@test Clapeyron.excess_gibbs_free_energy(system,p,T,z) ≈ Clapeyron.test_excess_gibbs_free_energy(system,p,T,z)  rtol = 1e-6
    end
    @testset "Wilson" begin
        system = Wilson(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.530046633499114 rtol = 1e-6
        @test Clapeyron.activity_coefficient(system,p,T,z) ≈ Clapeyron.test_activity_coefficient(system,p,T,z)  rtol = 1e-6
        @test Clapeyron.excess_gibbs_free_energy(system,p,T,z) ≈ Clapeyron.test_excess_gibbs_free_energy(system,p,T,z)  rtol = 1e-6
    end

    @testset "NRTL" begin
        system = NRTL(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5309354738922405 rtol = 1e-6
    end

    @testset "aspen-NRTL" begin
        nrtl_vanilla = NRTL(["methanol","benzene"])
        system = aspenNRTL(["methanol","benzene"])
        system2 = aspenNRTL(nrtl_vanilla)
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5309354738922405 rtol = 1e-6
        @test Clapeyron.activity_coefficient(system2,p,T,z)[1] ≈ 1.5309354738922405 rtol = 1e-6
    end

    @testset "UNIQUAC" begin
        system = UNIQUAC(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.3630421218486388 rtol = 1e-6
    end

    @testset "UNIFAC" begin
        system = UNIFAC(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5322232657797463 rtol = 1e-6
        #when fast UNIFAC works, it should pass this test.
        # system2 = UNIFAC(["methanol","benzene"])
        # prop2 = ()
        # @test Clapeyron.activity_coefficient(system2,1e-4,423.15,[0.,1.])  ≈ [2.0807335111878937,1.0] rtol = 1e-6
    end

    @testset "UNIFAC2" begin
        system = UNIFAC2([("acetaldehyde", ["CH3" => 1, "HCO" => 1]),("acetonitrile", ["CH3CN" => 1])]; puremodel=BasicIdeal())
        @test log(Clapeyron.activity_coefficient(system, NaN, 323.15, [0.5,0.5])[1]) ≈ 0.029527741236233 rtol = 1e-6
    end

    @testset "ogUNIFAC" begin
        system = ogUNIFAC(["methanol","benzene"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.5133696314734384 rtol = 1e-6
    end

    @testset "ogUNIFAC2" begin
        system = ogUNIFAC2([("R22",["HCCLF2" => 1]),("carbon disulfide",["CS2" => 1])], puremodel=BasicIdeal())
        @test log(Clapeyron.activity_coefficient(system, NaN, 298.15, [0.3693,0.6307])[1]) ≈ 0.613323250984226 rtol = 1e-6
    end

    @testset "UNIFAC-FV" begin
        system = UNIFACFV(["PMMA","PS"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 8.63962025235759 rtol = 1e-6
    end

    @testset "UNIFAC-FV-poly" begin
        system = UNIFACFVPoly(["PMMA","PS"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 4.8275769947121985 rtol = 1e-6
    end

    @testset "FH" begin
        system = FH(["PMMA","PS"],[100,100])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.8265799707907238 rtol = 1e-6
    end

    @testset "COSMOSAC02" begin
        system = COSMOSAC02(["water","ethanol"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.3871817962565904 rtol = 1e-6
        @test Clapeyron.excess_gibbs_free_energy(system,p,T,z) ≈ 610.5706657776052 rtol = 1e-6
    end

    @testset "COSMOSAC10" begin
        system = COSMOSAC10(["water","ethanol"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.4015660588643404 rtol = 1e-6
    end

    @testset "COSMOSACdsp" begin
        system = COSMOSACdsp(["water","ethanol"])
        @test Clapeyron.activity_coefficient(system,p,T,z)[1] ≈ 1.4398951117248127 rtol = 1e-6
    end
    end
    @printline
end

@testset "Ideal models" begin
    @printline
    let T = 298.15, V = 1e-4,p = 1e5,z = Clapeyron.SA[1.0],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "Joback" begin
        system = JobackIdeal(["hexane"])
        @test Clapeyron.VT_isobaric_heat_capacity(system,V,298.15) ≈ 143.22076150138616 rtol = 1e-6
        @test Clapeyron.crit_pure(system)[1] ≈ 500.2728274871347 rtol = 1e-6
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 9.210841420941021 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
        @test Clapeyron.mass_density(system,p,T,z) ≈ Clapeyron.molecular_weight(system,z)*p/(Rgas(system)*T)
        s0 = JobackIdeal("acetone")
        @test Clapeyron.JobackGC.T_c(s0)[1] ≈ 500.5590 rtol = 1e-6
        @test Clapeyron.JobackGC.P_c(s0)[1] ≈ 48.025e5 rtol = 1e-6
        @test Clapeyron.JobackGC.V_c(s0)[1] ≈ 209.5e-6 rtol = 1e-6
        @test Clapeyron.JobackGC.T_b(s0)[1] ≈ 322.1100 rtol = 1e-6
        @test Clapeyron.JobackGC.H_form(s0)[1] ≈ -217830.0 rtol = 1e-6
        @test Clapeyron.JobackGC.G_form(s0)[1] ≈ -154540.0 rtol = 1e-6
        @test Clapeyron.JobackGC.C_p(s0,300)[1] ≈ 75.3264 rtol = 1e-6
        @test Clapeyron.JobackGC.H_fusion(s0)[1] ≈ 5.1250e3 rtol = 1e-6
        @test Clapeyron.JobackGC.H_vap(s0)[1] ≈ 29.0180e3 rtol = 1e-6
        @test Clapeyron.JobackGC.Visc(s0,300)[1] ≈ 0.0002942 rtol = 9e-4
    end

    @testset "Reid" begin
        system = ReidIdeal(["butane"])
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 9.210842104089576 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
        @test Clapeyron.mass_density(system,p,T,z) ≈ Clapeyron.molecular_weight(system,z)*p/(Rgas(system)*T)
    end

    @testset "Shomate" begin
        system = ShomateIdeal(["water"])
        coeff = system.params.coeffs[1]
        @test Clapeyron.evalcoeff(system,coeff,500) ≈ 35.21836175 rtol = 1e-6
        @test Clapeyron.eval∫coeff(system,coeff,500) ≈ 15979.2447 rtol = 1e-6
        @test Clapeyron.eval∫coeffT(system,coeff,500) ≈ 191.00554 rtol = 1e-6
        @test Clapeyron.mass_density(system,p,T,z) ≈ Clapeyron.molecular_weight(system,z)*p/(Rgas(system)*T)
    end

    @testset "Walker" begin
        system = WalkerIdeal(["hexane"])
        @test Clapeyron.molecular_weight(system)*1000 ≈ 86.21
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 179.51502015696653 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
        @test Clapeyron.mass_density(system,p,T,z) ≈ Clapeyron.molecular_weight(system,z)*p/(Rgas(system)*T)
        
        ideal_csv = """
        Clapeyron Database File
        Walker Ideal Like Parameters  [csvtype = like,grouptype = Walker]
        species,Mw,Nrot,theta1,theta2,theta3,theta4,deg1,deg2,deg3,deg4,source
        ideal,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,ideal gas
        """
        system2 = WalkerIdeal(["ideal gas" => ["ideal" => 1]],userlocations = [ideal_csv])
        @test Clapeyron.isobaric_heat_capacity(system2,1e5,T)/Clapeyron.Rgas() ≈ 2.5 
    end

    @testset "Monomer" begin
        system = MonomerIdeal(["hexane"])
        @test Clapeyron.a_ideal(system,V,T,z) ≈ -10.00711774776317 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
        @test Clapeyron.mass_density(system,p,T,z) ≈ Clapeyron.molecular_weight(system,z)*p/(Rgas(system)*T)
    end

    @testset "Empiric" begin
        #Empiric Ideal from JSON
        system = EmpiricIdeal(["water"])
        #
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 7.932205569922042 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
        @test Clapeyron.mass_density(system,p,T,z) ≈ Clapeyron.molecular_weight(system,z)*p/(Rgas(system)*T)

        #Empiric Ideal from already existing MultiFluid model
        system = Clapeyron.idealmodel(MultiFluid(["water"]))
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 7.932205569922042 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14

        #Empiric Ideal from already existing single fluid model
        system = Clapeyron.idealmodel(system.pures[1])
        @test Clapeyron.a_ideal(system,V,T,z) ≈ 7.932205569922042 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
    end

    @testset "Aly-Lee" begin
        system = AlyLeeIdeal(["methane"])
        @test_broken Clapeyron.a_ideal(system,V,T,z) ≈ 9.239701647126086 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
        @test Clapeyron.mass_density(system,p,T,z) ≈ Clapeyron.molecular_weight(system,z)*p/(Rgas(system)*T)

        #we use the default GERG 2008 parameters for methane, test if the Cp is equal
        system_gerg = Clapeyron.idealmodel(GERG2008(["methane"]))
        Cp_system = Clapeyron.VT_isobaric_heat_capacity(system,V,T,z)
        Cp_gerg = Clapeyron.VT_isobaric_heat_capacity(system_gerg,V,T,z)

        @test Cp_system ≈ Cp_gerg rtol = 5e-5
    end

    @testset "Cp - LNG - Estimation" begin
        #Mw to obtain γ₀ = 0.708451
        system = CPLNGEstIdeal(["a1"],userlocations = (;Mw = [20.5200706797]))
        #test at 324.33 K, paper says Cp = 44.232, but the calculations in the paper seem off
        @test Clapeyron.VT_isobaric_heat_capacity(system,0.03,324.33) ≈ 44.231 rtol = 5e-4
    end

    @testset "PPDS" begin
        m1 = PPDSIdeal("krypton")
        @test isobaric_heat_capacity(m1,1,303.15)/Rgas(m1) ≈ 2.5
        @test Clapeyron.mass_density(m1,p,T,z) ≈ Clapeyron.molecular_weight(m1,z)*p/(Rgas(m1)*T)
        mw2 = 32.042 #MonomerIdeal("methanol").params.Mw.values[1]
        m2 = PPDSIdeal("methanol")
        #verification point in ref 1, table A.6
        @test isobaric_heat_capacity(m2,1,303.15)/mw2 ≈ 1.3840 rtol = 1e-4
    end
end
    @printline
end

@testset "Multi-parameter models" begin
    @printline
    let T = 298.15, V = 1e-4,p = 1e5,z = Clapeyron.SA[1.0],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333]; 
    #warning, we are in the pseudo maxwell loop, those properties are nonsense, but they evaluate anyway.
    @printline
    @testset "IAPWS95" begin
        system = IAPWS95()
        system_ideal = Clapeyron.idealmodel(system)
        @test Clapeyron.a_ideal(system_ideal, V, T, z) ≈ 7.9322055699220435 rtol = 1e-6
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 7.9322055699220435 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -2.1152889226862166e14 rtol = 1e-6
        #because we are in this regime, numerical accuracy suffers. that is why big(V) is used instead.
        @test Clapeyron.ideal_consistency(system,big(V),T,z) ≈ 0.0 atol = 1e-14
    end

    @testset "PropaneRef" begin
        system = PropaneRef()
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 0.6426994942361217 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -2.436280448227229 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z) ≈ 0.0 atol = 1e-14
    end

    @testset "GERG2008" begin
        system = GERG2008(["water"])
        @test Clapeyron.a_ideal(system, V, T, z) ≈ 4.500151936577565 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T, z) ≈ -10.122119572808764 rtol = 1e-6
        z4   = [0.25,0.25,0.25,0.25]
        system = GERG2008(["water","carbon dioxide","hydrogen sulfide","argon"])
        @test Clapeyron.a_ideal(system, V, T, z4) ≈ 3.1136322215343917 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system, V, T, z4) ≈ 0.0 rtol = 1e-14
        @test Clapeyron.a_res(system, V, T, z4) ≈ -1.1706377677539772 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,V,T,z4) ≈ 0.0 atol = 1e-14

        system_R = Clapeyron.GERG2008(["methane","ethane"],Rgas = 8.2)
        @test Clapeyron.Rgas(system_R) == 8.2
    end

    @testset "EOS-LNG" begin
        #LNG paper, table 16
        Tx = 150.0
        Vx = 1/(18002.169)
        zx   = [0.6,0.4]
        system = EOS_LNG(["methane","butane"])
        dep = departure_functions(system)
        @test count(!iszero,dep) == 1
        @test Clapeyron.eos(system,Vx,Tx,zx) ≈ -6020.0044 rtol = 5e-6

        #546
        system_546 = EOS_CG(["carbon dioxide","water"])
        @test volume(system_546,1e6,293.15,[0.1742127426126829, 0.8257872573873171]) ≈ 2.080669615020994e-5 rtol = 1e-6
    end

    @testset "LKP" begin
        Tx = 150.0
        Vx = 1/(18002.169)
        zx   = [0.6,0.4]
        system = LKP(["methane","butane"])
        test_scales(system)
        test_k(system)
        #Clapeyron.a_res(EOS_LNG(["methane","butane"]),V,T,z) ≈ -6.56838705236683
        @test Clapeyron.a_res(system,Vx,Tx,zx) ≈ -6.469596957611441 rtol = 5e-6
        system2 = LKPSJT(["methane","butane"])
        #TODO:check this value
        @test Clapeyron.a_res(system2,Vx,Tx,zx) ≈ -5.636923220762173 rtol = 5e-6
    end

    @testset "LJRef" begin
        system = LJRef(["methane"])
        Tx = 1.051*Clapeyron.T_scale(system)
        Vx = Clapeyron._v_scale(system)/0.673
        @test Clapeyron.a_ideal(system, Vx, Tx) ≈ 5.704213386278148 rtol = 1e-6
        @test Clapeyron.a_res(system, Vx, Tx) ≈ -2.244730279521925 rtol = 1e-6
        @test Clapeyron.ideal_consistency(system,Vx,Tx,z1) ≈ 0.0 atol = 1e-14
    end

    @testset "Xiang-Deiters" begin
        system = XiangDeiters(["water"])
        #equal to Clapeyron.a_ideal(BasicIdeal(["water"]), V, T, z)
        @test Clapeyron.a_ideal(system, V, T, z1) ≈ -0.33605470137749016 rtol = 1e-6
        @test Clapeyron.a_res(system, V, T)  ≈ -34.16747927719535 rtol = 1e-6
    end

    @testset "multiparameter misc" begin
        T = 300.0
        V = 1/200
        z2 = [0.5,0.5]
        z1 = Clapeyron.SA[1.0]
        #=
        apart from CO2+H2 (because we use a more recent departure), all models are compared with coolprop outputs:
        PropsSI("alphar","Dmolar|gas",200.0,"T",300.0,fluid)
        =#
        model = MultiFluid(["carbon dioxide","hydrogen"],verbose = true) #test verbose and gauss+exponential
        test_scales(model)
        test_repr(model.departure.params.parameters[1,2],str = ["Departure MultiParameter coefficients:","Fij:","Polynomial power terms:","Gaussian bell-shaped terms:"])
        test_repr(model.pures[1],str = ["Exponential terms: 27","Non Analytic terms:","Polynomial power terms:","Gaussian bell-shaped terms:"])

        pures = Clapeyron.split_pure_model(model)
        @test pures isa Vector{SingleFluid{EmpiricAncillary}}
        @test Clapeyron.wilson_k_values(model,1e6,300.0) ≈ [6.738566125478432, 54.26124873240438] rtol = 1e-3
        @test Clapeyron.a_res(model,V,T,z2) ≈ -0.005482930754339683 rtol = 1e-6
        model2 = SingleFluid("ammonia",verbose = true) #test Gaob parser
        @test Clapeyron.a_res(model2,V,T,z1) ≈ -0.05006143389915488 rtol = 1e-6
        model3 = SingleFluid("D4",verbose = true) #ideal CP0 parser
        @test Clapeyron.a_ideal(model3,V,T,z1) ≈ 0.8441669238992482 rtol = 1e-6
        model4 = SingleFluid("R14",verbose = true) #ResidualHelmholtzExponential
        @test Clapeyron.a_res(model4,V,T,z1) ≈ -0.017855323645451636 rtol = 1e-6
        model5 = SingleFluid("water",Rgas = 10.0)
        @test Rgas(model5) == 10.0
    end
    @printline
    end
end

@testset "SPUNG models" begin
    @printline
    let T = 298.15, V = 1e-4,p = 1e5,z = Clapeyron.SA[1.0],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333]; 
    @testset "SPUNG (SRK)" begin
        system = SPUNG(["ethane"])
        @test Clapeyron.shape_factors(system, V, T, z)[1] ≈ 0.824678830913322 rtol = 1e-6
    end


    @testset "SPUNG (PCSAFT)" begin
        system = SPUNG(["ethane"],PropaneRef(),PCSAFT(["ethane"]),PCSAFT(["propane"]))
        @test Clapeyron.shape_factors(system, V, T, z)[1] ≈ 0.8090183134644525 rtol = 1e-6
    end
    end
end

@testset "lattice models" begin
    @printline
    let T = 298.15, V = 1e-4,p = 1e5,z = Clapeyron.SA[1.0],z1 = Clapeyron.SA[1.0],z2 = [0.5,0.5],z3 = [0.333, 0.333,0.333]; 

    @testset "single component" begin
        system = Clapeyron.SanchezLacombe(["carbon dioxide"])
        @test Clapeyron.a_res(system, V, T, z) ≈ -0.9511044462267396 rtol = 1e-6
    end

    @testset "Sanchez-Lacombe,Kij rule" begin
        system = SanchezLacombe(["carbon dioxide","benzoic acid"],mixing = SLKRule)
        @test Clapeyron.a_res(system, V, T, z2) ≈ -6.494291842858994 rtol = 1e-6
    end

    @testset "Sanchez-Lacombe K0-K1-L rule" begin
        system = SanchezLacombe(["carbon dioxide","benzoic acid"],mixing = SLk0k1lMixingRule)
        @test Clapeyron.a_res(system, V, T, z2) ≈ -5.579621796375229 rtol = 1e-6
    end
    end
end


@testset "Correlations" begin
    @testset "Saturation" begin
        @testset "AntoineSat" begin
            system = AntoineEqSat(["water"])
            p0 = saturation_pressure(system,400.01)[1]
            @test p0 ≈ 244930.458389 rtol = 1e-6
            @test saturation_temperature(system,p0)[1] ≈ 400.01 rtol = 1e-6
        end

        @testset "DIPPR101Sat" begin
            system = DIPPR101Sat(["water"])
            p0 = saturation_pressure(system,400.01)[1]
            @test p0 ≈ 245338.15099198322 rtol = 1e-6
            @test saturation_temperature(system,p0)[1] ≈ 400.01 rtol = 1e-6
        end
        
        @testset "LeeKeslerSat" begin
            system = LeeKeslerSat(["water"])
            p0 = saturation_pressure(system,400.01)[1]
            @test p0 ≈ 231731.79240876858 rtol = 1e-6
            @test saturation_temperature(system,p0)[1] ≈ 400.01 rtol = 1e-6
        end

        #=@testset "PolExpSat" begin
            system = PolExpSat(["water"])
            p0 = saturation_pressure(system,400.01)[1]
            @test p0 ≈ 231731.79240876858 rtol = 1e-6
            @test saturation_temperature(system,p0)[1] ≈ 400.01 rtol = 1e-6
        end=#
    end
    @testset "LiquidVolume" begin
        @testset "COSTALD" begin
            system = COSTALD(["water"])
            @test volume(system,1e5,300.15) ≈ 1.8553472145724288e-5 rtol = 1e-6
            system2 = COSTALD(["water","methanol"])
            @test volume(system2,1e5,300.15,[0.5,0.5]) ≈ 2.834714146558056e-5 rtol = 1e-6
            @test volume(system2,1e5,300.15,[1.,0.]) ≈ 1.8553472145724288e-5 rtol = 1e-6
        end

        @testset "DIPPR105Liquid" begin
            system = DIPPR105Liquid(["water"])
            @test volume(system,1e5,300.15) ≈ 1.8057164858551208e-5 rtol = 1e-6
            system2 = DIPPR105Liquid(["water","methanol"])
            @test volume(system2,1e5,300.15,[0.5,0.5]) ≈ 2.9367802477146567e-5 rtol = 1e-6
            @test volume(system2,1e5,300.15,[1.,0.]) ≈ 1.8057164858551208e-5 rtol = 1e-6
        end

        @testset "RackettLiquid" begin
            system = RackettLiquid(["water"])
            @test volume(system,1e5,300.15) ≈ 1.6837207241594103e-5 rtol = 1e-6
            system2 = RackettLiquid(["water","methanol"])
            @test volume(system2,1e5,300.15,[0.5,0.5]) ≈ 3.2516352601748416e-5 rtol = 1e-6
            @test volume(system2,1e5,300.15,[1.,0.]) ≈ 1.6837207241594103e-5 rtol = 1e-6
        end

        @testset "YamadaGunnLiquid" begin
            system = YamadaGunnLiquid(["water"])
            @test volume(system,1e5,300.15) ≈ 2.0757546189420953e-5 rtol = 1e-6
            system2 = YamadaGunnLiquid(["water","methanol"])
            @test volume(system2,1e5,300.15,[0.5,0.5]) ≈ 3.825994563177142e-5 rtol = 1e-6
            @test volume(system2,1e5,300.15,[1.,0.]) ≈ 2.0757546189420953e-5 rtol = 1e-6
        end

        @testset "Grenke-Elliott water" begin
            system = GrenkeElliottWater()
            system_ref = IAPWS95()
            p1,T1 = 611.657, 273.16
            p2,T2 = 101325.0, 273.152519
            for Ti in range(250.0,280.0,30)
                for log10Pi in range(5,8,20)
                    Pi = exp10(log10Pi)
                    v = volume(system,Pi,Ti)
                    v_ref = volume(system_ref,Pi,Ti,phase = :l)
                    @test v ≈ v_ref rtol = 1e-3
                    @test pressure(system,v,Ti) ≈ pressure(system_ref,v_ref,Ti) rtol = 1e-3
                end
                Cpi = Clapeyron.mass_isobaric_heat_capacity(system,101325.0,Ti)
                Cpi_test = Clapeyron.water_cp(system,Ti)
                @test Cpi ≈ Cpi_test rtol = 1e-6
            end
        end

        @testset "Holten Water" begin
            model = HoltenWater()
            Tc = 228.2
            Rm = 461.523087
            ρ0 = 1081.6482

            #table 8
            verification_data = [
                273.15 0.101325 999.84229  −0.683042  5.088499  4218.3002 1402.3886 0.09665472 0.62120474
                235.15 0.101325 968.09999  −29.633816 11.580785 5997.5632 1134.5855 0.25510286 0.091763676
                250    200      1090.45677 3.267768   3.361311  3708.3902 1668.2020 0.03042927 0.72377081
                200    400      1185.02800 6.716009   2.567237  3338.525  1899.3294 0.00717008 1.1553965
                250    400      1151.71517 4.929927   2.277029  3757.2144 2015.8782 0.00535884 1.4345145
                ]
            
            for i in 1:5
                T = verification_data[i,1]
                p = verification_data[i,2]*1e6
                @test mass_density(model,p,T) ≈ verification_data[i,3] rtol = 1e-6
                @test isobaric_expansivity(model,p,T)*1e4 ≈ verification_data[i,4] rtol = 1e-6
                @test isothermal_compressibility(model,p,T)*1e10 ≈ verification_data[i,5] rtol = 1e-6
                @test mass_isobaric_heat_capacity(model,p,T) ≈ verification_data[i,6] rtol = 1e-6
                @test speed_of_sound(model,p,T) ≈ verification_data[i,7] rtol = 1e-6
                

                𝕡 = p/(Rm*Tc*ρ0)
                L = Clapeyron.water_L(model,𝕡,T)
                ω = 2 + 0.5212269*𝕡
                xe = Clapeyron.water_x_frac(model,L,ω)
                @test xe ≈ verification_data[i,8] rtol = 1e-6
                @test L ≈ verification_data[i,9] rtol = 1e-6
            end
        end
    end

    @testset "Virial Coeff" begin
        @testset "AbbottVirial" begin
            system = AbbottVirial(["methane","ethane"])
            @test volume(system,1e5,300,[0.5,0.5]) ≈ 0.024820060368027988 rtol = 1e-6
            #a_res(PR,0.05,300,[0.5,0.5]) == -0.0023705490820905483
            @test Clapeyron.a_res(system,0.05,300,[0.5,0.5]) ≈ -0.0024543543773067693 rtol = 1e-6
        end

        @testset "TsonopoulosVirial" begin
            system = TsonopoulosVirial(["methane","ethane"])
            @test volume(system,1e5,300,[0.5,0.5]) ≈ 0.02485310667780686 rtol = 1e-6
            #a_res(PR,0.05,300,[0.5,0.5]) == -0.0023705490820905483
            @test Clapeyron.a_res(system,0.05,300,[0.5,0.5]) ≈ -0.0017990881811592349 rtol = 1e-6
        end

        @testset "EoSVirial2" begin
            cub = PR(["methane","ethane"])
            system = EoSVirial2(cub)
            #exact equality here, as cubics have an exact second virial coefficient
            @test volume(system,1e5,300,[0.5,0.5]) == Clapeyron.volume_virial(cub,1e5,300,[0.5,0.5])
            #a_res(PR,0.05,300,[0.5,0.5]) == -0.0023705490820905483
            @test Clapeyron.a_res(system,0.05,300,[0.5,0.5]) ≈ -0.0023728381262076137 rtol = 1e-6
        end
    end

    @testset "Solid Models" begin
        @testset "SolidHfus" begin
            model = SolidHfus(["water"])
            @test chemical_potential(model,1e5,298.15,[1.])[1] ≈ 549.1488193300384 rtol = 1e-6
        end

        @testset "SolidKs" begin
            model = SolidKs(["water"])
            @test chemical_potential(model,1e5,298.15,[1.])[1] ≈ 549.1488193300384 rtol = 1e-6
        end

        @testset "IAPWS-06" begin
            model = IAPWS06()
            #table 6 of Ice-Rev2009 document
            p1,T1 = 611.657, 273.16
            p2,T2 = 101325.0, 273.152519
            p3,T3 = 100e6, 100.0
            Mw = Clapeyron.molecular_weight(model)
            @test mass_gibbs_energy(model,p1,T1) ≈ 0.611784135 rtol = 1e-6
            @test mass_gibbs_energy(model,p2,T2) ≈ 0.10134274069e3 rtol = 1e-6
            @test mass_gibbs_energy(model,p3,T3) ≈ -0.222296513088e6 rtol = 1e-6
            @test volume(model,p1,T1)/Mw ≈ 0.109085812737e-2 rtol = 1e-6
            @test volume(model,p2,T2)/Mw ≈ 0.109084388214e-2 rtol = 1e-6
            @test volume(model,p3,T3)/Mw ≈ 0.106193389260e-2 rtol = 1e-6
            test_volume(model,p1,T1)
            test_volume(model,p2,T2)
            test_volume(model,p3,T3)
            @test mass_isobaric_heat_capacity(model,p1,T1) ≈ 0.209678431622e4 rtol = 1e-6
            @test mass_isobaric_heat_capacity(model,p2,T2) ≈ 0.209671391024e4 rtol = 1e-6
            @test mass_isobaric_heat_capacity(model,p3,T3) ≈ 0.866333195517e3 rtol = 1e-6
            @test isentropic_compressibility(model,p1,T1) ≈ 0.114161597779e-9 rtol = 1e-6
            @test isentropic_compressibility(model,p2,T2) ≈ 0.114154442556e-9 rtol = 1e-6
            @test isentropic_compressibility(model,p3,T3) ≈ 0.886060982687e-10 rtol = 1e-6
            @test isothermal_compressibility(model,p1,T1) ≈ 0.117793449348e-9 rtol = 1e-6
            @test isothermal_compressibility(model,p2,T2) ≈ 0.117785291765e-9 rtol = 1e-6
            @test isothermal_compressibility(model,p3,T3) ≈ 0.886880048115e-10 rtol = 1e-6
        end

        @testset "Jäger-Span solid CO2" begin
            model = JagerSpanSolidCO2()
            #table 5 of 10.1021/je2011677 
            p1,T1 = 0.51795e6, 216.592
            p2,T2 = 100e6, 100.0
            @test gibbs_energy(model,p1,T1) ≈ -1.447007522e3 rtol = 1e-6
            @test gibbs_energy(model,p2,T2) ≈ -2.961795962e3 rtol = 1e-6
            @test volume(model,p1,T1) ≈ 2.848595255e-5 rtol = 1e-6
            @test volume(model,p2,T2) ≈ 2.614596591e-5 rtol = 1e-6
            test_volume(model,p1,T1)
            test_volume(model,p2,T2)
            @test entropy(model,p1,T1) ≈ -1.803247012e1 rtol = 1e-6
            @test entropy(model,p2,T2) ≈ -5.623154438e1 rtol = 1e-6
            @test isobaric_heat_capacity(model,p1,T1) ≈ 5.913420271e1 rtol = 1e-6
            @test isobaric_heat_capacity(model,p2,T2) ≈ 3.911045710e1 rtol = 1e-6
            @test isobaric_expansivity(model,p1,T1) ≈ 8.127788321e-4 rtol = 1e-6
            @test isobaric_expansivity(model,p2,T2) ≈ 3.843376525e-4 rtol = 1e-6
            @test isothermal_compressibility(model,p1,T1) ≈ 2.813585169e-10 rtol = 1e-6
            @test isothermal_compressibility(model,p2,T2) ≈ 1.149061787e-10 rtol = 1e-6
        end
    end
end
