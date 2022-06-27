@testset "utils" begin
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
    @testset "split_model" begin
        models2 = split_model(model2)
        @test models2[1].components[1] == model2.components[1]
        @test models2[2].components[1] == model2.components[2]
        @test models2[1].icomponents == models2[2].icomponents == 1:1

        model2_unsplit = only(split_model(model2,[[1,2]]))
        @test model2_unsplit.icomponents == model2.icomponents
        @test model2_unsplit.components == model2.components

        model4_split = Clapeyron.split_model(model4)
        @test model4_split[4].groups.n_groups[1][2] == 3
        @test model4_split[1].components[1] == "methane"
        @test all(isone(length(model4_split[i].components)) for i in 1:4)

        gc3_split = Clapeyron.split_model(gc3)
        @test all(isone(length(gc3_split[i].components)) for i in 1:3)
        @test all(isone(length(gc3_split[i].puremodel)) for i in 1:3)
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
        @test repr("text/plain",ideal1) == "WalkerIdeal{BasicIdeal} with 1 component:\n \"hexane\": \"CH3\" => 2, \"CH2\" => 4\nContains parameters: Mw, Nrot, theta1, theta2, theta3, theta4, deg1, deg2, deg3, deg4"
        #@newmodel
        @test repr(model2) == "PCSAFT{BasicIdeal}(\"water\", \"ethanol\")"
        @test repr("text/plain",model2) == "PCSAFT{BasicIdeal} with 2 components:\n \"water\"\n \"ethanol\"\nContains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol"
        #@newmodelsimple
        @test repr(noparam1) == "NoTranslation()"
        @test repr("text/plain",noparam1) == "NoTranslation\n"
        @test repr(simple1) == "PRAlpha(\"propane\")"
        @test repr("text/plain",simple1) == "PRAlpha with 1 component:\n \"propane\"\nContains parameters: acentricfactor"
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
    end
    @printline

end
