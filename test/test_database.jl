using OpenSAFT, Test

filepath_like = "test_csvs/data_test_like.csv"
filepath_unlike = "test_csvs/data_test_unlike.csv"
filepath_assoc = "test_csvs/data_test_assoc.csv"

@testset "csv_parser" begin
    @testset "parseline and tryparse_fallback" begin
        parsedline = OpenSAFT.parseline(filepath_like, 4)
        @test typeof(parsedline) <: Array
        @test parsedline[1] == "sp1"
        @test typeof(parsedline[2]) == Missing
        @test typeof(parsedline[3]) == Float64
    end
    @testset "findmatches" begin
        found_lines = OpenSAFT.findmatches(filepath_like, "16", 3)
        @test found_lines == OpenSAFT.findmatches(filepath_like, "16", "Mw"; header_row = 3)
        @test found_lines == [4, 9]
        found_lines = OpenSAFT.findmatches(filepath_like, "nonexisting_entry", 3)
        @test found_lines == []
        @test_throws ErrorException OpenSAFT.findmatches(filepath_like, "16", "nonexisting_header"; header_row = 3)
    end
    @testset "findmatches_pair" begin
        found_lines = OpenSAFT.findmatches_pair(filepath_assoc, "sp4", "sp5", 1, 4)
        @test found_lines == OpenSAFT.findmatches_pair(filepath_assoc, "sp4", "sp5", "species1", "species2"; header_row=3)
        @test found_lines == [9, 10]
        @test_throws ErrorException OpenSAFT.findmatches_pair(filepath_assoc, "sp4", "sp5", "species1", "nonexisting_header"; header_row = 1)
    end
end

@testset "extractdatabase" begin
    @testset "auxilliary" begin
        @testset "models_return" begin
            models = OpenSAFT.models_return("None")
            @test typeof(models) == Array{String, 1}
            @test "PCSAFT" in models
            models = OpenSAFT.models_return("PCSAFT")
            @test models == ["PCSAFT"]
            @test_throws ErrorException OpenSAFT.models_return("nonexstant_model")
        end
        @testset "createfilepath" begin
            filepath = OpenSAFT.createfilepath("PCSAFT", "like") # use a variant when supported
            @test filepath == dirname(pathof(OpenSAFT)) * "/../database/PCSAFT/data_PCSAFT_like.csv"
        end
        @testset "customdatabase_check" begin
            @test_throws ErrorException OpenSAFT.customdatabase_check("None", "like", filepath_like)
        end
    end
    @testset "searchdatabase" begin
        @testset "searchdatabase_like" begin
            found_models = OpenSAFT.searchdatabase_like(["ethanol", "water", "nonexisting_species"])
            @test haskey(found_models[Set(["water"])], "PCSAFT")
            @test haskey(found_models, Set(["nonexisting_species"]))
            @test isempty(found_models[Set(["nonexisting_species"])])

            found_models = OpenSAFT.searchdatabase_like(["ethanol", "water"], "PCSAFT")
            @test haskey(found_models[Set(["water"])], "PCSAFT")
            @test_throws ErrorException OpenSAFT.searchdatabase_like(["ethanol", "nonexisting_species"], "PCSAFT")

            found_models = OpenSAFT.searchdatabase_like(["sp1"], "PCSAFT"; customdatabase = filepath_like) 
            @test found_models[Set(["sp1"])]["PCSAFT"] == 9
            @test_throws ErrorException OpenSAFT.searchdatabase_like(["sp1", "nonexistant_species"], "PCSAFT"; customdatabase = filepath_like) 
        end
        @testset "searchdatabase_unlike" begin
            found_models = OpenSAFT.searchdatabase_unlike(["methane", "butane", "nonexisting_species"])
            @test haskey(found_models[Set(["methane", "butane"])], "PCSAFT")
            @test !haskey(found_models, Set(["methane", "nonexistant_species"]))

            found_models = OpenSAFT.searchdatabase_unlike(["methane", "butane", "nonexistant_species"], "PCSAFT")
            @test haskey(found_models[Set(["methane", "butane"])], "PCSAFT")
            @test !haskey(found_models, Set(["methane", "nonexistant_species"]))

            found_models = OpenSAFT.searchdatabase_unlike(["sp1", "sp2", "nonexistant_species"], "PCSAFT"; customdatabase = filepath_unlike) 
            @test found_models[Set(["sp1", "sp2"])]["PCSAFT"] == 7
            @test !haskey(found_models, Set(["sp1", "nonexistant_species"]))
        end
        @testset "searchdatabase_assoc" begin
            found_models = OpenSAFT.searchdatabase_assoc(["methanol", "ethanol", "nonexistant_species"])
            @test haskey(found_models[Set(["methanol"])], "PCSAFT")
            @test !haskey(found_models, Set(["nonexistant_species"]))

            found_models = OpenSAFT.searchdatabase_assoc(["methanol", "ethanol", "nonexistant_species"], "PCSAFT")
            @test haskey(found_models[Set(["methanol"])], "PCSAFT")
            @test !haskey(found_models, Set(["methanol", "nonexistant_species"]))

            found_models = OpenSAFT.searchdatabase_assoc(["sp4", "sp5", "nonexistant_species"], "PCSAFT"; customdatabase = filepath_assoc) 
            @test found_models[Set(["sp4", "sp5"])]["PCSAFT"] == [9, 10]
            @test !haskey(found_models, Set(["sp1", "nonexistant_species"]))
        end
    end
    @testset "retrieveparams" begin
        @testset "retrieveparams_like" begin
            found_params = OpenSAFT.retrieveparams_like(["ethanol", "water"], "PCSAFT")
            @test found_params[Set(["ethanol"])]["sigma"] == 3.1771
        end
        @testset "retrieveparams_unlike" begin
            found_params = OpenSAFT.retrieveparams_unlike(["methane", "butane", "pentane"], "PCSAFT")
            @test found_params[Set(["methane", "butane"])]["k"] == 0.022
            @test found_params[Set(["pentane", "methane"])]["k"] == 0.024
        end
        @testset "retrieveparams_assoc" begin
            found_params = OpenSAFT.retrieveparams_assoc(["methanol", "ethanol"], "PCSAFT")
            @test found_params[Set([(Set(["methanol"]), "e"), (Set(["methanol"]), "H")])]["epsilon_assoc"] == 2899.5
            found_params = OpenSAFT.retrieveparams_assoc(["sp2", "sp3", "sp4", "sp5"], "PCSAFT"; customdatabase = filepath_assoc)
            @test found_params[Set([(Set(["sp3"]), "e"), (Set(["sp3"]), "H")])]["epsilon_assoc"] == 2700
            @test found_params[Set([(Set(["sp5"]), "e"), (Set(["sp5"]), "H")])]["epsilon_assoc"] == 2100
            @test found_params[Set([(Set(["sp5"]), "e")])]["epsilon_assoc"] == 2200
            @test found_params[Set([(Set(["sp3"]), "H"), (Set(["sp4"]), "H")])]["epsilon_assoc"] == 2300
            @test found_params[Set([(Set(["sp3"]), "H"), (Set(["sp4"]), "e")])]["epsilon_assoc"] == 2400
            @test !haskey(found_params, Set([(Set(["sp3"]), "e"), (Set(["sp4"]), "H")]))
            @test found_params[Set([(Set(["sp4"]), "H"), (Set(["sp5"]), "e")])]["epsilon_assoc"] == 2500
            @test found_params[Set([(Set(["sp4"]), "e"), (Set(["sp5"]), "H")])]["epsilon_assoc"] == 2600
        end
    end
end

@testset "filterparams" begin
    raw_params = OpenSAFT.retrieveparams(["sp1", "sp2", "sp3", "sp4", "sp5"], "PCSAFT";
        customdatabase_like=filepath_like, customdatabase_unlike=filepath_unlike, customdatabase_assoc=filepath_assoc)
    like_params_dict, unlike_params_dict, assoc_params_dict =
        OpenSAFT.filterparams(raw_params, ["epsilon", "m"]; unlike_params = ["k"], assoc_params = ["epsilon_assoc", "bond_vol"])
    @test like_params_dict["m"][Set(["sp1"])] == 1.0
    @test unlike_params_dict["k"][Set(["sp1", "sp2"])] == 0.02
    @test assoc_params_dict["epsilon_assoc"][Set([(Set(["sp3"]), "H"), (Set(["sp3"]), "e")])] == 2700
    @test_throws KeyError OpenSAFT.filterparams(raw_params, ["nonexisting_parameter"])
    @test_throws KeyError OpenSAFT.filterparams(raw_params, ["m"]; assoc_params = ["nonexisting_parameter"])
end
