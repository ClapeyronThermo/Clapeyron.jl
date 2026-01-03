@testset "Reported errors" begin
    #https://github.com/ClapeyronThermo/Clapeyron.jl/issues/104
    @testset "#104" begin
        model = VTPR(["carbon dioxide"])
        p = 1e5
        T = 273.15
        @test fugacity_coefficient(model, p, T)[1] ≈ 0.9928755806091429 rtol = 1E-6
        @test activity_coefficient(model, p, T)[1] ≈ 1.0
    end
    @testset "#112" begin
        model = CPA(["methanol"])
        @test crit_pure(model)[1] ≈ 538.2329369300235 rtol = 1e-6
    end

    @testset "DM - SAFTgammaMie 1" begin
        #this constructor was failing on Clapeyron, 3.10-dev
        model=SAFTgammaMie(["water","ethyl acetate"])
        assocparam = model.vrmodel.params.bondvol
        @test assocparam.sites[1] == ["H2O/H", "H2O/e1"]
        @test assocparam.sites[2] == ["COO/e1"]


        #a symmetric matrix was generated from non-symmetric GC values, Clapeyron 0.6.4
        model_mix = SAFTgammaMie([("MEA",["NH2"=>1]),("Carbon Dioxide",["CO2"=>1])];
            userlocations = (Mw = [16.02285, 44.01],
            epsilon = [284.78 134.58;
                        134.58 207.89],
            sigma = [3.2477, 3.05],
            lambda_a = [6, 5.055],
            lambda_r = [10.354 50.06;
                        50.06  26.408],
            vst = [1, 2],
            S = [0.79675, 0.84680],
            n_H=[2, 0],
            n_e=[1, 0],
            n_a1=[0, 1],
            n_a2=[0, 1],
            epsilon_assoc = Dict([(("NH2","H"),("NH2","e")) => 1070.80,
                                    (("CO2","a1"),("NH2","e")) => 3313,
                                    (("CO2","a2"),("NH2","e")) => 4943.6]),
            bondvol = Dict([(("NH2","H"),("NH2","e")) => 95.225e-30,
                            (("CO2","a1"),("NH2","e")) => 3280.3e-30,
                            (("CO2","a2"),("NH2","e")) => 142.64e-30])))


        bondvol_mixed = model_mix.vrmodel.params.bondvol
        co2 = "Carbon Dioxide"
        mea = "MEA"
        #normal
        @test bondvol_mixed[(co2,"CO2/a1"),(mea,"NH2/e")] == 3280.3e-30
        @test bondvol_mixed[(co2,"CO2/a2"),(mea,"NH2/e")] == 142.64e-30
        #reverse
        @test bondvol_mixed[(mea,"NH2/e"),(co2,"CO2/a1")] == 3280.3e-30
        @test bondvol_mixed[(mea,"NH2/e"),(co2,"CO2/a2")] == 142.64e-30
        
        #swich sites: this should result in zeros:
        @test bondvol_mixed[(co2,"CO2/e"),(mea,"NH2/a1")] == 0
        @test bondvol_mixed[(mea,"CO2/e"),(co2,"NH2/a2")] == 0
        @test bondvol_mixed[(mea,"NH2/a1"),(co2,"CO2/e")] == 0
        @test bondvol_mixed[(co2,"NH2/a2"),(mea,"CO2/e")] == 0
    end

    @testset "#140" begin
    #FractionVector with length(x) == 0.
        model = PCSAFT(["water","carbon dioxide"])
        res = bubble_pressure(model,280,Clapeyron.FractionVector(0.01),ChemPotBubblePressure(nonvolatiles = ["water"]))
        @test res[1] ≈ 4.0772545187410433e6 rtol = 1e-6
    end

    @testset "#145" begin
        #incorrect indexing of gc_to_comp_sites when there is more than one assoc site.
        like_data = """
        Clapeyron Database File
        SAFTgammaMie Like Parameters [csvtype = like,grouptype = SAFTgammaMie]
        species,vst,S,lambda_r,lambda_a,sigma,epsilon,n_H,n_e1,n_e2,Mw
        CH2_PEO,1,1,12,6,4,300,0,0,0,100
        cO_1sit,1,1,12,6,4,300,0,1,0,100
        """

        mw_data = """
        Clapeyron Database File
        SAFTgammaMie Like Parameters
        species,Mw
        CH2_PEO,100
        cO_1sit,100
        """

        assoc_data = """
        Clapeyron Database File,,,,,,
        SAFTgammaMie Assoc Parameters [csvtype = assoc,grouptype = SAFTgammaMie]
        species1,site1,species2,site2,epsilon_assoc,bondvol,source
        H2O,H,cO_1sit,e1,2193.2,5e-29,
        CH2OH,H,cO_1sit,e1,1572.5,4.331e-28,
        """

        group_data = """
        Clapeyron Database File,
        SAFTgammaMie Groups [csvtype = groups,grouptype = SAFTgammaMie]
        species,groups
        PEG_1sit,"[""CH2_PEO"" => 2000,""cO_1sit"" => 1000,""CH2OH"" => 2]"
        """

        model = SAFTgammaMie(["water","PEG_1sit"],userlocations = [like_data,assoc_data,mw_data],group_userlocations = [group_data])
        #in this case,because the groups are 1 to 1 correspondence to each molecule, the amount of groups should be the same
        @test length(model.vrmodel.params.epsilon_assoc.values.values) == length(model.params.epsilon_assoc.values.values)
        #test if we got the number of sites right
        @test model.vrmodel.sites.n_sites[2][1] == 1000 #1000 sites cO_1sit/e1 in PEG.
    end

    @testset "#154" begin
        #there was a problem when using the @newmodel macros outside the Clapeyron module. this should suffice as a test.
        abstract type PCSAFTModel_test <: SAFTModel end

        # Defining the parameters used by the model
        struct PCSAFTParam_test <: EoSParam
            Mw::SingleParam{Float64}
            segment::SingleParam{Float64}
            sigma::PairParam{Float64}
            epsilon::PairParam{Float64}
            epsilon_assoc::AssocParam{Float64}
            bondvol::AssocParam{Float64}
        end

        # Creating a model struct called PCSAFT, which is a sub-type of PCSAFTModel, and uses parameters defined in PCSAFTParam
        @newmodel PCSAFT_test PCSAFTModel_test PCSAFTParam_test
        @newmodelsimple PCSAFT_testsimple PCSAFTModel_test PCSAFTParam_test
        @newmodelgc PCSAFT_testgc PCSAFTModel_test PCSAFTParam_test

        #test if the macros actually generate something
        @test PCSAFT_test <: EoSModel #@newmodel
        @test PCSAFT_testsimple <: EoSModel #@newmodelsimple
        @test PCSAFT_testgc <: EoSModel #@newmodelgc
    end


    @testset "#162" begin
        #a longstanding problem, init_model didn't worked with functions.
        #a long time ago, SRK was a model, but now it is just a function that returns an RK model.
        model1 = Wilson(["water","ethanol"];puremodel=SRK)
        @test model1 isa Clapeyron.EoSModel

        #this case is just for compatibility with the notebooks that were originally released.
        model2 = VTPR(["carbon monoxide","carbon dioxide"];alpha = BMAlpha)
        @test model2 isa Clapeyron.EoSModel
    end

    @testset "#171, #366" begin
        #=
        #171
        This is a problem that occurs in an intersection between split_model and cross-association sites
        a single component model created from scratch don't have any cross association sites,
        but a single component model created from split_model does have those sites.
        we neet to check if the results of both are equal.
        =#
        model = SAFTgammaMie(["water","acetone"])
        model_split = split_model(model)[2]
        model_pure = SAFTgammaMie(["acetone"])
        res_pure = Clapeyron.eos(model_pure,1.013e6,298.15) #works
        res_split = Clapeyron.eos(model_split,1.013e6,298.15) #should work
        @test res_pure ≈ res_split

        #=
        #366
        incorrect conversion of MixedGCSegmentParam.
        =#
        mix_segment_f64 = model.params.mixed_segment
        mix_segment_bigfloat = convert(Clapeyron.MixedGCSegmentParam{BigFloat},mix_segment_f64)
        @test mix_segment_f64.values.v ≈ mix_segment_bigfloat.values.v
    end

    @testset "#188" begin
        #=
        in cubics, we do a separation between the alpha model and the main model. while this allows for unprecedented
        flexibility, this also complicates the case of the default model. there is also an unrelated error, about needing to pass Vc,
        because our database has all the critical parameters in one file
        =#
        data = (
            species = ["A", "B"],
            Tc = [18.0, 3.5],
            Pc = [2.3, 0.1],
            Mw = [1.0, 1.0],
            acentricfactor = [0.1, 0.3]
        )

        file = ParamTable(
            :single,
            data,
            name="db1"
        )

        system = PR(["A", "B"], userlocations = [file])
        @test system.params.Tc[2] == 3.5

    end

    @testset "#201" begin
        #=
        Ternary LLE
        =#
        if hasfield(aspenNRTL,:puremodel)
            @test aspenNRTL(["water", "acetone", "dichloromethane"],puremodel = PR) isa EoSModel
        else
            @test aspenNRTL(["water", "acetone", "dichloromethane"]) isa EoSModel
        end

        @test UNIFAC(["water", "acetone", "dichloromethane"]) isa EoSModel
    end

    @testset "#212" begin
        #=
        Polar PCSAFT
        it uses `a` and `b` as site names
        =#
        model = PCPSAFT("water")
        @test model isa EoSModel
        @test "a" in model.sites.flattenedsites
        @test "b" in model.sites.flattenedsites
    end

    @testset "#75 - regression" begin
        #=
        pharmaPCSAFT, incorrect mixing rules
        =#

        userlocations = (Mw = [352.77,18.01],
        segment = [14.174,1.2046817736],
        sigma = [3.372,2.7927],
        epsilon = [221.261,353.9449],
        n_H = [2,1],
        n_e = [2,1],
        k = [0 -0.155;-0.155 0],
        kT = [0 2.641e-4;2.641e-4 0],
        epsilon_assoc = Dict(
        (("griseofulvin","e"),("griseofulvin","H")) => 1985.49,
        (("water08","e"),("water08","H")) => 2425.6714),
        bondvol = Dict(
        (("griseofulvin","e"),("griseofulvin","H")) => 0.02,
        (("water08","e"),("water08","H")) => 0.04509)
        )
        model = pharmaPCSAFT(["griseofulvin","water08"];userlocations = userlocations)
        #found while testing, test that we don't mutate the input
        @test userlocations.sigma == [3.372,2.7927]
        T = 310.15
        p = 1.01325e5
        z=[7.54e-7, 1-7.54e-7]

        γ1 = activity_coefficient(model,p,T,z)

        #this was an error too. check commit that added this
        #@test γ1[1] ≈ 55334.605821130834 rtol = 1e-4

        @test γ1[1] ≈ 51930.06908022231 rtol = 1e-4

    end

    @testset "#262" begin
        #=
        SAFT-gamma Mie
        Ensure ij,ab == ji,ba, and that the group-comp is correct for assoc
        =#
        model = SAFTgammaMie(["water","ethanol"])
        @test model.params.epsilon_assoc.values[1,2][1,2] == model.params.epsilon_assoc.values[2,1][2,1]
        @test model.params.epsilon_assoc.values[1,2][1,2] == model.vrmodel.params.epsilon_assoc.values[1,2][1,2]
    end

    @testset "SorptionModels.jl - init kij with user" begin
        #=
        on SL, passing k in userlocations did not work.
        =#


        v★(P★, T★,) = 8.31446261815324 * T★ / P★ / 1000000 # J / (mol*K) * K / mpa -> pa * m3 / (mol * mpa) ->  need to divide by 1000000 to get m3/mol
        ϵ★(T★) = 8.31446261815324 * T★ # J / (mol * K) * K -> J/mol
        r(P★, T★, ρ★, mw) = mw * (P★ * 1000000) / (8.31446261815324 * T★ * (ρ★ / 1e-6)) # g/mol * mpa * 1000000 pa/mpa / ((j/mol*K) * K * g/(cm3 / 1e-6 m3/cm3)) -> unitless

        P★ = [534., 630.]
        T★ = [755., 300.]
        ρ★ = [1.275, 1.515]
        mw = [100000, 44.01]
        kij = [0 -0.0005; -0.0005 0]
        model1 = Clapeyron.SL(
            ["PC", "CO2"],
            userlocations = Dict(
                :vol => v★.(P★, T★,),
                :segment => r.(P★, T★, ρ★, mw),
                :epsilon => ϵ★.(T★),
                :Mw => mw
            ),
            mixing_userlocations = (;k0 = kij, k1 = [0 0; 0 0], l = [0 0; 0 0])
        )

        model2 = Clapeyron.SL(
            ["PC", "CO2"],
            userlocations = Dict(
                :vol => v★.(P★, T★,),
                :segment => r.(P★, T★, ρ★, mw),
                :epsilon => ϵ★.(T★),
                :Mw => mw,
                :k => kij
            )
        )
        @test Clapeyron.get_k(model1)[1] ≈ Clapeyron.get_k(model2)[1]
    end

    @testset "https://github.com/ClapeyronThermo/Clapeyron.jl/discussions/239" begin
        #test for easier initialization of CPA/SAFT without association
        m1 = Clapeyron.CPA(["Methanol"])
        m2 = CPA(["Methanol"]; userlocations=(;
        a = m1.params.a.values[1],
        b = m1.params.b.values[1],
        c1 = m1.params.c1.values,
        Mw = m1.params.Mw.values,
        Tc = m1.params.Tc.values,
        Pc = m1.cubicmodel.params.Pc.values,
        n_H = [1],
        n_e = [1],
        epsilon_assoc = Dict((("Methanol","H"),("Methanol","e")) => m1.params.epsilon_assoc.values.values[1]),
        bondvol = Dict((("Methanol","H"),("Methanol","e")) => m1.params.bondvol.values.values[1]))
        )
        @test volume(m1,1e5,333.0) ≈ volume(m2,1e5,333.0)

        m3 = PCSAFT("water",userlocations =(segment = 1,Mw = 1,epsilon = 1,sigma = 1.0))
        @test length(m3.params.bondvol.values.values) == 0
    end

    @testset "infinite dilution derivatives - association" begin
        #reported by viviviena in slack
        cp_params = (a = [36.54206320678348, 39.19437678197436, 25.7415, 34.91747774761048], b = [-0.03480434051958945, -0.05808483585041852, 0.2355, -0.014935581577635826], c = [0.000116818199785053, 0.0003501220208504329, 0.0001578, 0.000756101594841365], d = [-1.3003819534791665e-7, -3.6941157412454843e-7, -4.0939e-7, -1.0894144551347726e-6], e = [5.2547403746728466e-11, 1.276270011886522e-10, 2.1166e-10, 4.896983427747592e-10])
        idealmodel = ReidIdeal(["water", "methanol", "propyleneglycol","methyloxirane"]; userlocations = cp_params)
        pcpsaft = PCPSAFT(["water", "methanol", "propyleneglycol","methyloxirane"], idealmodel = idealmodel)
        water = split_model(pcpsaft)[1]
        z = [1.0,0.0,0.0,0.0]
        v = volume(pcpsaft, 101325.0, 298.15, z, phase = :liquid)
        cp_pure = Clapeyron.VT_isobaric_heat_capacity(water, v, 298.15, 1.0)
        cp_mix = Clapeyron.VT_isobaric_heat_capacity(pcpsaft, v, 298.15, z)
        @test cp_mix ≈ cp_pure
        @test cp_mix ≈ 69.21259493306137
    end

    @testset "CPA init" begin
        model1 = CPA(["C3","C4"];
        userlocations=(;
        Mw = [44.097,58.124],
        Tc = [369.83,425.18],
        a = [9.11875,13.14274],
        b = [0.057834,0.072081],
        c1 = [0.6307,0.70771],
        epsilon_assoc=nothing,
        bondvol=nothing))
        @test model1 isa CPA
    end

    @testset "#357 - electrolyte equilibria" begin

        model1 = ePCSAFT(["water"], ["calcium", "chloride"])
        salts1 = [("calcium chloride", ("calcium" => 1, "chloride" => 2))]
        x1 = molality_to_composition(model1, salts1, 1.0)
        bub_test = 374.7581484748338
        bubP1_test = 2971.744917038001

        bubT1 = bubble_temperature(model1, 101325, x1, FugBubbleTemperature(nonvolatiles = ["calcium", "chloride"]))[1]
        @test bubT1 ≈ bub_test rtol = 1e-6

        bubT1_chempot = bubble_temperature(model1, 101325, x1, ChemPotBubbleTemperature(T0=373.0, nonvolatiles=["calcium","chloride"]))[1]
        @test_broken bubT1_chempot ≈ bub_test rtol = 1e-6

        bubT1_chempot2 = bubble_temperature(model1, 101325, x1, ChemPotBubbleTemperature(nonvolatiles=["calcium","chloride"]))[1]
        @test_broken bubT1_chempot2 ≈ bub_test rtol = 1e-6

        bubT1_2 = bubble_temperature(model1, 101325, x1, FugBubbleTemperature(nonvolatiles = ["calcium", "chloride"],y0=[1.,0.,0.],vol0=(1.8e-5,1.),T0=373.15))[1]
        @test bubT1_2 ≈ bub_test rtol = 1e-6

        bubP1 = bubble_pressure(model1,298.15,x1,FugBubblePressure(nonvolatiles=["calcium","chloride"],y0=[1.,0.,0.],vol0=(1e-5,1.)))[1]
        @test bubP1 ≈ bubP1_test rtol = 1e-6

        bubP1_2 = bubble_pressure(model1,298.15,x1,FugBubblePressure(nonvolatiles=["calcium","chloride"]))[1]
        @test bubP1 ≈ bubP1_test rtol = 1e-6
    end

    @testset "#416" begin
        #if we build a cubic and only provide critical parameters (without acentric factor), read the database
        #to build the alpha model
        model = PR(["nitrogen"]; userlocations=(;
           Tc = [0.3],
           Pc = [1.],
           Mw = [1.]))
        @test !model.alpha.params.acentricfactor.ismissingvalues[1]
    end

    @testset "#513" begin
        #=
        error in x0_sat_pure_spinodal
        incorrect bounds, causing failure in convergence.
        =#
        model = tcRK("R1243zf")
        Tr = range(0.2,1.0,500)
        Tc = first(crit_pure(model))
        psat = first.(saturation_pressure.(model,Tr .* Tc))
        @test count(isnan,psat[100:end]) == 0
    end
end
