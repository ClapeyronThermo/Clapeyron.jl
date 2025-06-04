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
        @info "The following 2 error messages are expected:"
        @test_throws ArgumentError split_model(noparam1)
        @test_throws MethodError split_model(model2,missing)
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

        noparam1_split = split_model(noparam1,1:5)
        @test length(noparam1_split) == 5
        @test noparam1_split[1] == noparam1

        #splitting parameters
        @test split_model(model2.params.segment)[1][1] == model2.params.segment[1]
        @test split_model(model2.params.sigma)[1][1,1] == model2.params.sigma[1,1]

        #from notebooks, #173
        nb_test = SAFTgammaMie(["methane","nitrogen","carbon dioxide","ethane","propane","butane","isobutane",
        "pentane","isopentane","hexane","heptane","octane"])
        @test length(split_model(nb_test)) == 12

        #weird error found on splitting groups
        model0 = SAFTgammaMie(["ethane"])
        model0_split = SAFTgammaMie(["methane","ethane"]) |> split_model |> last
        @test model0.params.epsilon.values[1,1] == model0_split.params.epsilon.values[1,1]
    end

    @testset "export_model" begin
        @testset "SAFT Model" begin
            model_og = PCSAFT(["water","ethanol"])
            export_model(model_og)
            model_ex = PCSAFT(["water","ethanol"]; userlocations = ["singledata_PCSAFT.csv","pairdata_PCSAFT.csv","assocdata_PCSAFT.csv"])

            @test model_og.params.segment.values == model_ex.params.segment.values
            @test model_og.params.epsilon.values == model_ex.params.epsilon.values
            @test model_og.params.epsilon_assoc.values.values == model_ex.params.epsilon_assoc.values.values
        end

        @testset "Cubic Model" begin
            model_og = PR(["water","ethanol"])
            export_model(model_og)
            model_ex = PR(["water","ethanol"]; userlocations = ["singledata_PR.csv","pairdata_PR.csv"],
                                       alpha_userlocations = ["singledata_PRAlpha.csv"])

            @test model_og.params.a.values == model_ex.params.a.values
            @test model_og.alpha.params.acentricfactor.values == model_ex.alpha.params.acentricfactor.values
        end

        @testset "Activity & GC Model" begin
            model_og = UNIFAC(["water","ethanol"])
            export_model(model_og)
            model_ex = UNIFAC(["water","ethanol"]; userlocations = ["singledata_UNIFAC.csv","pairdata_UNIFAC.csv"])

            @test model_og.params.Q.values == model_ex.params.Q.values
            @test model_og.params.A.values == model_ex.params.A.values
        end

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

    
    @testset "has_sites-has_groups" begin
        @test has_sites(typeof(gc3)) == false
        @test has_sites(typeof(model2)) == has_sites(typeof(model4)) == has_sites(typeof(gc2)) == true
        @test has_groups(typeof(gc2)) == has_groups(typeof(gc3)) == has_groups(typeof(ideal1)) == true
        @test has_groups(typeof(simple1)) == has_groups(typeof(model2)) == false
    end

    @testset "eosshow" begin
        #@newmodelgc
        @test repr(ideal1) == "WalkerIdeal{BasicIdeal}(\"hexane\")"
        @test repr("text/plain",ideal1) == "WalkerIdeal{BasicIdeal} with 1 component:\n \"hexane\": \"CH3\" => 2, \"CH2\" => 4\nGroup Type: Walker\nContains parameters: Mw, Nrot, theta1, theta2, theta3, theta4, deg1, deg2, deg3, deg4, reference_state"
        #@newmodel
        @test repr(model2) == "PCSAFT{BasicIdeal, Float64}(\"water\", \"ethanol\")"
        @test repr("text/plain",model2) == "PCSAFT{BasicIdeal, Float64} with 2 components:\n \"water\"\n \"ethanol\"\nContains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol"
        #@newmodelsimple
        @test repr(noparam1) == "NoTranslation()"
        @test repr("text/plain",noparam1) == "NoTranslation()"
        @test repr(simple1) == "PRAlpha(\"propane\")"
        @test repr("text/plain",simple1) == "PRAlpha with 1 component:\n \"propane\"\nContains parameters: acentricfactor"
    end

    @testset "Clapeyron Param show" begin
        @test repr(model2.params) == "Clapeyron.PCSAFTParam{Float64}"
        @test repr("text/plain",model2.params) == "Clapeyron.PCSAFTParam{Float64} for [\"water\", \"ethanol\"] with 6 params:\n Mw::SingleParam{Float64}\n segment::SingleParam{Float64}\n sigma::PairParam{Float64}\n epsilon::PairParam{Float64}\n epsilon_assoc::AssocParam{Float64}\n bondvol::AssocParam{Float64}"
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
        @test startswith(Clapeyron.doi2bib("10.1021/I160057A011"),"@article{Peng_1976")

    end

    @testset "alternative input" begin
        @test PCSAFT("water" => ["H2O"=>1],idealmodel = WalkerIdeal) isa EoSModel
        @test PCSAFT(["water" => ["H2O"=>1]],idealmodel = WalkerIdeal) isa EoSModel
        @test PCSAFT("water") isa EoSModel
        @test JobackIdeal("hexane") isa EoSModel 
        @test PCSAFT(["water" => ["H2O"=>1]]) isa EoSModel
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

    @printline
    if Base.VERSION >= v"1.8" #for some reason, it segfaults on julia 1.6
        @testset "ambiguities" begin
            ambiguities = Test.detect_ambiguities(Clapeyron)
            if !isempty(ambiguities)
                foreach(display, ambiguities)
            end
            @test length(ambiguities) == 0
        end
    end
 end
