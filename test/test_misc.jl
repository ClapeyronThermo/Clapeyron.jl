#struct for testset "#161"
struct PCSAFT161 <: Clapeyron.PCSAFTModel
    components::Vector{String}
    params::Clapeyron.PCSAFTParam
    references::Vector{String}
    weird_thing::Int
end

@testset "misc" begin
    @printline
    model2 = PCSAFT(["water","ethanol"])
    model4 = SAFTgammaMie(["methane","butane","isobutane","pentane"])
    gc3 = UNIFAC(["propane","butane","isobutane"])
    gc2 = SAFTgammaMie([
        "ethanol",
        ("ibuprofen", ["CH3"=>3, "COOH"=>1, "aCCH"=>1, "aCCH2"=>1, "aCH"=>4])])

    ideal1 = WalkerIdeal(["hexane"])
    noparam1 = gc3.puremodel[1].translation
    simple1 = gc3.puremodel[1].alpha
    model_structgc = structSAFTgammaMie(["ethanol","octane"])
    @testset "split_model" begin
        models2 = split_model(model2)
        @test models2[1].components[1] == model2.components[1]
        @test models2[2].components[1] == model2.components[2]

        model2_unsplit = only(split_model(model2,[[1,2]]))
        @test model2_unsplit.components == model2.components

        model4_split = Clapeyron.split_model(model4)
        @test model4_split[4].groups.n_groups[1][2] == 3
        @test model4_split[1].components[1] == "methane"
        @test all(isone(length(model4_split[i])) for i in 1:4)

        gc3_split = Clapeyron.split_model(gc3)
        @test all(isone(length(gc3_split[i])) for i in 1:3)
        @test all(isone(length(gc3_split[i].puremodel)) for i in 1:3)

        structgc_split = Clapeyron.split_model(model_structgc)
        @test structgc_split[1].groups.n_intergroups[1] == [0 1; 1 0]
        @test structgc_split[2].groups.n_intergroups[1] == [0 2; 2 5]
    end

    @testset "single component error" begin
        model = PCSAFT(["water","methane"])
        @test_throws DimensionMismatch saturation_pressure(model,300.15)
        @test_throws DimensionMismatch crit_pure(model)
        @test_throws DimensionMismatch saturation_temperature(model,1e5)
        @test_throws DimensionMismatch acentric_factor(model)
        @test_throws DimensionMismatch enthalpy_vap(model,300.15)
        @test_throws DimensionMismatch Clapeyron.x0_sat_pure(model,300.15)
        @test_throws DimensionMismatch saturation_liquid_density(model,300.15)
    end
    
    @testset "macros" begin
        comps(model) = Clapeyron.@comps
        groups(model) = Clapeyron.@groups
        groups(model,i) = Clapeyron.@groups(i)
        sites(model,i) = Clapeyron.@sites(i)
        f_eos(model,V,T,z,c) = 2 + c
        model = V = T = z = nothing

        @test groups(gc2) == 1:6
        @test comps(gc2) == 1:2
        @test groups(gc2,2) == [2,3,4,5,6]
        @test groups(gc3,1) == groups(gc3,2) #propane and butane has the same amount of groups
        @test sites(gc2.vrmodel,2) == [3,4,5]
        @test sites(model2,1) == [1,2]
        @test Clapeyron.@f(f_eos,pi) == 2+pi
        @test Clapeyron.@nan(Base.log(-1),3) == 3
        @test_throws MethodError Clapeyron.@nan(Base.log("s"),3)
    end

    using Clapeyron: has_sites,has_groups
    @testset "has_sites-has_groups" begin
        @test has_sites(typeof(gc3)) == false
        @test has_sites(typeof(model2)) == has_sites(typeof(model4)) == has_sites(typeof(gc2)) == true
        @test has_groups(typeof(gc2)) == has_groups(typeof(gc3)) == has_groups(typeof(ideal1)) == true
        @test has_groups(typeof(simple1)) == has_groups(typeof(model2)) == false
    end

    @testset "eosshow" begin
        #@newmodelgc
        @test repr(ideal1) == "WalkerIdeal{BasicIdeal}(\"hexane\")"
        @test repr("text/plain",ideal1) == "WalkerIdeal{BasicIdeal} with 1 component:\n \"hexane\": \"CH3\" => 2, \"CH2\" => 4\nGroup Type: Walker\nContains parameters: Mw, Nrot, theta1, theta2, theta3, theta4, deg1, deg2, deg3, deg4"
        #@newmodel
        @test repr(model2) == "PCSAFT{BasicIdeal}(\"water\", \"ethanol\")"
        @test repr("text/plain",model2) == "PCSAFT{BasicIdeal} with 2 components:\n \"water\"\n \"ethanol\"\nContains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol"
        #@newmodelsimple
        @test repr(noparam1) == "NoTranslation()"
        @test repr("text/plain",noparam1) == "NoTranslation\n"
        @test repr(simple1) == "PRAlpha(\"propane\")"
        @test repr("text/plain",simple1) == "PRAlpha with 1 component:\n \"propane\"\nContains parameters: acentricfactor"
    end

    @testset "Clapeyron Param show" begin
        @test repr(model2.params) == "Clapeyron.PCSAFTParam"
        @test repr("text/plain",model2.params) == "Clapeyron.PCSAFTParam for [\"water\", \"ethanol\"] with 6 params:\n Mw::SingleParam{Float64}\n segment::SingleParam{Float64}\n sigma::PairParam{Float64}\n epsilon::PairParam{Float64}\n epsilon_assoc::AssocParam{Float64}\n bondvol::AssocParam{Float64}"
    end

    @testset "phase symbols" begin
        @test Clapeyron.canonical_phase(:l) == :liquid
        @test Clapeyron.canonical_phase(:v) == :vapour
        @test Clapeyron.canonical_phase(:m) == :m
    end

    @testset "citing" begin
        umr = UMRPR(["water"],idealmodel = WalkerIdeal)
        #citations = ["10.1021/I160057A011", "10.1021/ie049580p", "10.1021/i260064a004", "10.1021/acs.jced.0c00723"] |> Set
        citation_full = Clapeyron.cite(umr) |> Set
        citation_top = Clapeyron.doi(umr) |> Set
        citation_ideal = Clapeyron.cite(umr.idealmodel) |> Set
        citation_mixing = Clapeyron.cite(umr.mixing) |> Set
        citation_translation = Clapeyron.cite(umr.translation) |> Set
        @test citation_ideal ⊆ citation_full
        @test citation_top ⊆ citation_full
        @test citation_mixing ⊆ citation_full
        @test citation_translation ⊆ citation_full
        _io = Base.IOBuffer()
        Clapeyron.show_references(_io,umr)
        citation_show = String(take!(_io))
        @test citation_show == "\nReferences: 10.1021/I160057A011, 10.1021/ie049580p, 10.1021/i260064a004, 10.1021/acs.jced.0c00723"
        @test Clapeyron.doi2bib("10.1021/I160057A011") == "@article{Peng_1976,\n\tdoi = {10.1021/i160057a011},\n\turl = {https://doi.org/10.1021%2Fi160057a011},\n\tyear = 1976,\n\tmonth = {feb},\n\tpublisher = {American Chemical Society ({ACS})},\n\tvolume = {15},\n\tnumber = {1},\n\tpages = {59--64},\n\tauthor = {Ding-Yu Peng and Donald B. Robinson},\n\ttitle = {A New Two-Constant Equation of State},\n\tjournal = {Industrial {\\&}amp\$\\mathsemicolon\$ Engineering Chemistry Fundamentals}\n}"
    end
    @printline

    @testset "core utils" begin
        @test Clapeyron.parameterless_type(typeof(rand(5))) === Array
        @test Clapeyron._vecparser("1 2 3") == [1,2,3]
        @test Clapeyron._vecparser("1 2 3.5") == [1,2,3.5]
        @test_throws ErrorException Clapeyron._vecparser("not numbers")
        @test Clapeyron.split_2("a b") == ("a","b")
        @test Clapeyron.split_2("a|b",'|') == ("a","b")
    end

    @testset "Reported errors" begin
        #https://github.com/ypaul21/Clapeyron.jl/issues/104
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

    end
    @printline
    if Base.VERSION >= v"1.8" #for some reason, it segfaults on julia 1.6
        @testset "ambiguities" begin
            ambiguities = Test.detect_ambiguities(Clapeyron)
            @test length(ambiguities) == 0
        end
    end
    #testset for equilibria bugs
    @testset "challenging equilibria" begin
        #@testset "dew_temperature N°1" begin
        #    modelp = PCSAFT(["water","methanol"])
        #end
    end
 end
