@testset "Reported errors" begin
    #https://github.com/ClapeyronThermo/Clapeyron.jl/issues/104
    @testset "#104" begin
        model = VTPR(["carbon dioxide"])
        p = 1e5
        T = 273.15
        @test fugacity_coefficient(model, p, T)[1] ≈ 0.9928244080356565 rtol = 1E-6
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

    @testset "#161" begin
        #problems with registermodel
        Clapeyron.@registermodel PCSAFT161
        @test hasmethod(Base.length,Tuple{PCSAFT161})
        @test hasmethod(Base.show,Tuple{IO,PCSAFT161})
        @test hasmethod(Base.show,Tuple{IO,MIME"text/plain",PCSAFT161})
        @test hasmethod(Clapeyron.molecular_weight,Tuple{PCSAFT161,Array{Float64}})
    end

    @testset "#162" begin
        #a longstanding problem, init_model didn't worked with functions.
        #a long time ago, SRK was a model, but now it is just a function that returns an RK model.
        model1 = Wilson(["water","ethanol"];puremodel=SRK)
        @test model1 isa Clapeyron.EoSModel

        #this case is just for compatibility with the notebooks that were originally released.
        model2 = VTPR(["carbon monoxide","carbon dioxide"];alpha=BMAlpha)
        @test model2 isa Clapeyron.EoSModel
    end

    @testset "#171" begin
        #=
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
        @test aspenNRTL(["water", "acetone", "dichloromethane"],puremodel = PR) isa EoSModel
        @test UNIFAC(["water", "acetone", "dichloromethane"]) isa EoSModel

    end

    @testset "#212" begin
        #=
        Polar PCSAFT
        it uses `a` and `b` as site names
        =#
        model = PPCSAFT("water")
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
        
        T = 310.15
        p = 1.01325e5
        z=[7.54e-7, 1-7.54e-7]

        γ1 = activity_coefficient(model,p,T,z)

        @test γ1[1] ≈ 55334.605821130834 rtol = 1e-4
    end
end
