using Clapeyron, Test, LinearAlgebra

@testset "Database Parsing" begin

    @testset "normalisestring" begin
        #before, it was coo, making a collision on Electrolyte SAFTgammaMie
        @test Clapeyron.normalisestring("COO-") == "coo-"
    end
     
    # The rest of the test will be conducted with a custom dataset in the test_csvs directory.
    testspecies = ["sp1", "sp2", "sp3", "sp4", "sp5"]

    filepath_normal = ["test_csvs/normal_pair_test.csv",
                        "test_csvs/normal_single_test.csv",
                        "test_csvs/normal_assoc_test.csv",
                        "test_csvs/normal_single2_test.csv",
                        "test_csvs/normal_assoc2_test.csv"]

    filepath_clashingheaders = ["test_csvs/headercollision_single_test.csv",
                                "test_csvs/headercollision_assoc_test.csv"]

    filepath_asymmetry = ["test_csvs/asymmetry_pair_test.csv",
                            "test_csvs/asymmetry_assoc_test.csv"]

    filepath_multiple_identifiers = ["test_csvs/multiple_identifiers_single_test.csv",
                                        "test_csvs/multiple_identifiers_pair_test.csv",
                                        "test_csvs/multiple_identifiers_assoc_test.csv"]
    filepath_gc = ["test_csvs/group_test.csv"]
    filepath_param_gc = ["test_csvs/group_param_test.csv"]
    # Check that it detects the right sites.

    opts = Clapeyron.ParamOptions()
    opts2 = Clapeyron.ParamOptions(ignore_missing_singleparams=["emptyparam","missingparam"])
    allparams,allnotfoundparams = Clapeyron.createparams(testspecies, filepath_normal, opts) #merge all found params
    
    @testset "getparams - general" begin

        params1 = Clapeyron.getparams(["water", "methanol"], ["SAFT/PCSAFT"],return_sites=false)
        #test that we have a sigma inside params1
        @test haskey(params1, "sigma")
        
        for param_i in values(params1)
            test_repr(param_i)
        end
        #this fails, because there are missing params
        @test_throws MissingException Clapeyron.compile_params(testspecies,allparams,allnotfoundparams,nothing,opts) #generate ClapeyronParams
        
        # Check that it throws an error if ignore_missing_singleparams is not set to true.
        @test_throws MissingException Clapeyron.getparams(testspecies; userlocations=filepath_normal,return_sites = false)
    
        # Clashing headers between association and non-association parameters are not allowed
        @test_throws ErrorException Clapeyron.getparams(testspecies; userlocations=filepath_clashingheaders)

        # If parameter is not tagged as ignore_missing_singleparams, incomplete diagonal will throw an error
        #ignore_missing_singleparams=["asymmetricpair"]
        @test_throws MissingException Clapeyron.getparams(testspecies; userlocations=filepath_asymmetry,return_sites = false)


        # Also, since a non-missing value exists on the diagonal of "asymmetricpair",
        # and the diagonal contains missing values, it should throw an error
        # when parameter is not in ignore_missing_singleparams
        @test_throws MissingException Clapeyron.getparams(testspecies; userlocations=filepath_asymmetry, asymmetricparams=["asymmetricpair", "asymmetricassoc"])

    end
    
    @printline
    sites2 = Clapeyron.buildsites(testspecies,allparams,allnotfoundparams,opts2)
    result = Clapeyron.compile_params(testspecies,allparams,allnotfoundparams,sites2,opts2) #generate ClapeyronParams
    @testset "params - sites" begin
        @test Set.(sites2.sites) == [Set{String}(),
                                    Set{String}(),
                                    Set(["e", "e2", "H"]),
                                    Set(["e", "H"]),
                                    Set(["e", "e2", "H"])]

    end
    params = Clapeyron.getparams(testspecies; userlocations=filepath_normal, ignore_missing_singleparams=["emptyparam","missingparam"], verbose = true)
    sites = Clapeyron.SiteParam(params["intparam"].components)
    
    @testset "params - printing" begin
        
        #Printing: SingleParam
        @test repr(params["intparam"]) == "SingleParam{Int64}(\"intparam\")[\"sp1\", \"sp2\", \"sp3\", \"sp4\", \"sp5\"]"
        @test repr("text/plain",params["intparam"]) == "SingleParam{Int64}(\"intparam\") with 5 components:\n \"sp1\" => 6\n \"sp2\" => 2\n \"sp3\" => 7\n \"sp4\" => 4\n \"sp5\" => 5"
        
        
        #Printing: PairParam
        @test repr(params["overwriteparam"]) == "PairParam{Float64}(\"overwriteparam\")[\"sp1\", \"sp2\", \"sp3\", \"sp4\", \"sp5\"]"
        @test repr("text/plain",params["overwriteparam"]) == "5×5 PairParam{Float64}([\"sp1\", \"sp2\", \"sp3\", \"sp4\", \"sp5\"]) with values:\n 1.6  4.0  0.0  0.0  0.0\n 4.0  1.2  0.0  3.0  0.0\n 0.0  0.0  1.3  2.0  0.0\n 0.0  3.0  2.0  1.4  0.0\n 0.0  0.0  0.0  0.0  1.5"
        #Printing: AssocParam
        #TODO: fix
        #@test repr(params["overwriteassocparam"]) == "AssocParam{String}(\"overwriteassocparam\")[\"val1\", \"val8\", \"val5\", \"val4\", \"val7\", \"val6\", \"val3\", \"42\"]"
        #@test repr("text/plain",params["overwriteassocparam"]) == "AssocParam{String}([\"sp1\", \"sp2\", \"sp3\", \"sp4\", \"sp5\"]) with 8 values:\n(\"sp3\", \"e\") >=< (\"sp3\", \"H\"): val1\n(\"sp3\", \"e2\") >=< (\"sp3\", \"H\"): val8\n(\"sp3\", \"H\") >=< (\"sp4\", \"e\"): val5\n(\"sp3\", \"H\") >=< (\"sp4\", \"H\"): val4\n(\"sp4\", \"e\") >=< (\"sp5\", \"H\"): val7\n(\"sp4\", \"H\") >=< (\"sp5\", \"e2\"): val6\n(\"sp5\", \"e\") >=< (\"sp5\", \"e\"): val3\n(\"sp5\", \"e\") >=< (\"sp5\", \"H\"): 42"
        #Printing: SiteParam
        @test repr(sites) == "SiteParam[\"sp1\" => [], \"sp2\" => [], \"sp3\" => [], \"sp4\" => [], \"sp5\" => []]"
        @test repr("text/plain",sites) == "SiteParam with 5 components:\n \"sp1\": (no sites)\n \"sp2\": (no sites)\n \"sp3\": (no sites)\n \"sp4\": (no sites)\n \"sp5\": (no sites)"
        # Check that all the types are correct.
    end

    @testset "params - types and values" begin
        
        @test typeof(params["intparam"]) <: Clapeyron.SingleParam{Int}
        @test typeof(params["doubleparam"]) <: Clapeyron.SingleParam{Float64}
        @test typeof(params["stringparam"]) <: Clapeyron.SingleParam{String}
        
        # If column has both strings and numbers, they should all be strings.
        @test typeof(params["mixedparam"]) <: Clapeyron.SingleParam{String}

        # Contains missing values
        @test typeof(params["missingparam"]) <: Clapeyron.SingleParam{Int}

        # All missing values
        #before returned Clapeyron.SingleParam{Any}
        #now returns Clapeyron.SingleParam{Float64}
        @test typeof(params["emptyparam"]) <: Clapeyron.SingleParam{Float64}

        # Overwrite Int with Float, and also an upgrade from single to pair param
        @test typeof(params["overwriteparam"]) <: Clapeyron.PairParam{Float64}
        #test for missingness in PairParam, i got bitten by this before
        @test params["overwriteparam"].ismissingvalues ==  Bool[0 0 1 1 1;
                                                                0 0 1 0 1;
                                                                1 1 0 0 1;
                                                                1 0 0 0 1;
                                                                1 1 1 1 0]

        #@test typeof(params["boolparam"]) <: Clapeyron.SingleParam{Bool}
        #we can now convert directly, no need to parse diferently

        # Overwrite String with Int
        @test typeof(params["overwritestringparam"]) <: Clapeyron.SingleParam{String}
        # Overwrite Int with String
        @test typeof(params["overwriteassocparam"]) <: Clapeyron.AssocParam{String}
        # Check that values of "sp1" and "sp3" has been correctly overwritten.
        # "sp1" was overwritten from the same file
        # "sp3" was overwritten from a separate file, "normal_single2_test.csv"
        @test params["intparam"].values == [6, 2, 7, 4, 5]
        # Also check that "testsource1" and "testsource3" are not present, and in their place
        # we have "testsource6" and "testsource19".
        @test "testsource1" ∉ params["intparam"].sources
        @test "testsource3" ∉ params["intparam"].sources
        @test "testsource5" ∈ params["intparam"].sources
        @test "testsource19" ∈ params["intparam"].sources

        # Check that missing values have been correctly defaulted.
        #sp1 appears twice. one with one value, another with missing.
        #the idea is not to overwrite the existing value with a a missing
        @test params["missingparam"].values == [1, 2, 3, 0, 0]
        @test params["missingparam"].ismissingvalues == Bool[0, 0, 0, 1, 1]

        #test that Passing a single param to a pair param doesnt erase the missings
        singletopairparam = PairParam(params["missingparam"],"singletopairparam")
        diagmissingvalues = view(singletopairparam.ismissingvalues,1:6:25)
        @test all(diagmissingvalues .== params["missingparam"].ismissingvalues)
        diagmissingvalues .= true
        @test all(singletopairparam.ismissingvalues)

        # Check that the promotion form 1D to 2D array is succesful, with non-diagonal values present and symmetrical.
        @test params["overwriteparam"].values == [1.6  4.0  0.0  0.0  0.0
                                                4.0  1.2  0.0  3.0  0.0
                                                0.0  0.0  1.3  2.0  0.0
                                                0.0  3.0  2.0  1.4  0.0
                                                0.0  0.0  0.0  0.0  1.5]

        @test Clapeyron.diagvalues(params["overwriteparam"]) == [1.6, 1.2, 1.3, 1.4, 1.5]
        @test Clapeyron.diagvalues(1.23) == 1.23
         
        assoc_param = params["assocparam"]
        #diagonal 3-3
        @test assoc_param[("sp3","H"),("sp3","e")] == 2000
        @test assoc_param[("sp3","e"),("sp3","H")] == 2000
        @test assoc_param[("sp3","H"),("sp3","e2")] == 2700
        @test assoc_param[("sp3","e2"),("sp3","H")] == 2700

        #asym: 3-4
        @test assoc_param[("sp3","H"),("sp4","e")] == 2400
        @test assoc_param[("sp4","e"),("sp3","H")] == 2400
        @test assoc_param[("sp3","H"),("sp4","H")] == 2300
        @test assoc_param[("sp4","H"),("sp3","H")] == 2300
        @test assoc_param[("sp3","e"),("sp4","H")] == 0
        @test assoc_param[("sp4","H"),("sp3","e")] == 0
        @test assoc_param[("sp3","e"),("sp4","e")] == 0
        @test assoc_param[("sp4","e"),("sp3","e")] == 0
        @test size(assoc_param["sp3","sp4"]) == size(transpose(assoc_param["sp4","sp3"]))

        #asym 4-5
        @test assoc_param[("sp4","H"),("sp5","e2")] == 2500
        @test assoc_param[("sp5","e2"),("sp4","H")] == 2500
        @test assoc_param[("sp4","e2"),("sp5","H")] == 0
        @test assoc_param[("sp5","H"),("sp4","e2")] == 0

        #asym 4-5
        @test assoc_param[("sp4","H"),("sp5","e")] == 0
        @test assoc_param[("sp5","e"),("sp4","H")] == 0
        @test assoc_param[("sp4","e"),("sp5","H")] == 2600
        @test assoc_param[("sp5","H"),("sp4","e")] == 2600

        #diagonal 5-5
        @test assoc_param[("sp5","H"),("sp5","e")] == 2100
        @test assoc_param[("sp5","e"),("sp5","H")] == 2100
        @test assoc_param[("sp5","e"),("sp5","e")] == 2200

        #=
        assoc_param_values =
        [[Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,3)          ]  [Array{Int64}(undef,0,2)]  [Array{Int64}(undef,0,3)       ]
        [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,3)          ]  [Array{Int64}(undef,0,2)]  [Array{Int64}(undef,0,3)       ]
        [Array{Int64}(undef,3,0)]  [Array{Int64}(undef,3,0)]  [[0 0 2000; 0 0 2700; 2000 2700 0]]  [[0 0; 0 0; 2400 2300]  ]  [[0 0 0; 0 0 0; 0 0 0]         ]
        [Array{Int64}(undef,2,0)]  [Array{Int64}(undef,2,0)]  [[0 0 2400; 0 0 2300]             ]  [[0 0; 0 0]             ]  [[0 0 2600; 0 2500 0]          ]
        [Array{Int64}(undef,3,0)]  [Array{Int64}(undef,3,0)]  [[0 0 0; 0 0 0; 0 0 0]            ]  [[0 0; 0 2500; 2600 0]  ]  [[2200 0 2100; 0 0 0; 2100 0 0]]]
        =#
        asymmetricparams = Clapeyron.getparams(testspecies; userlocations=filepath_asymmetry, asymmetricparams=["asymmetricpair", "asymmetricassoc"], ignore_missing_singleparams=["asymmetricpair"])

        @test asymmetricparams["asymmetricpair"].values == [0.06  0.04  0.0  0.0   0.0
                                                      0.05  0.0   0.0  0.0   0.0
                                                      0.0   0.0   0.0  0.02  0.0
                                                      0.0   0.03  0.0  0.0   0.0
                                                      0.0   0.0   0.0  0.0   0.0]
        #=
        asymmetricparams["asymmetricassoc"].values ==
        [[Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,3)       ]  [Array{Int64}(undef,0,2) ]  [Array{Int64}(undef,0,3)    ]
         [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,3)       ]  [Array{Int64}(undef,0,2) ]  [Array{Int64}(undef,0,3)    ]
         [Array{Int64}(undef,3,0)]  [Array{Int64}(undef,3,0)]  [[0 0 2800; 0 0 0; 2000 2700 0]]  [[0 2900; 0 0; 2400 2300]]  [[0 0 0; 0 0 0; 0 0 0]      ]
         [Array{Int64}(undef,2,0)]  [Array{Int64}(undef,2,0)]  [[0 0 0; 0 0 0]                ]  [[0 0; 0 0]              ]  [[0 0 0; 0 2500 0]          ]
         [Array{Int64}(undef,3,0)]  [Array{Int64}(undef,3,0)]  [[0 0 0; 0 0 0; 0 0 0]         ]  [[0 0; 0 0; 2600 0]      ]  [[2200 0 0; 0 0 0; 2100 0 0]]]
        =#
        # Testing for multiple identifiers
        multiple_identifiers = Clapeyron.getparams(testspecies; userlocations=filepath_multiple_identifiers)
        @test multiple_identifiers["param_single"].values == [100, 200, 200, 300, 300]
        @test multiple_identifiers["param_pair"].values   == [1000  4000     0     0     0
                                                             4000  2000     0  6000     0
                                                                0     0  2000  5000     0
                                                                0  6000  5000  3000     0
                                                                0     0     0     0  3000]
        @test multiple_identifiers["param_assoc"].values[3,2][1,1] == 10000 &&
                multiple_identifiers["param_assoc"].values[5,1][1,1] == 20000    
    
    end


    @testset "params - conversion,broadcasting,indexing" begin
        #single param conversion and broadcasting
        comps = ["aa","bb"]
        floatbool = SingleParam("float to bool",comps,[1.0,0.0])
        intbool = SingleParam("int to bool",comps,[1,0])
        str_single = SingleParam("int to bool",comps,["111","222"])
        substr_single = SingleParam("int to bool",comps,chop.(["111","222"]))

        @test size(floatbool) == (2,)
        #SingleParam - conversion
        @test convert(SingleParam{Bool},floatbool) isa SingleParam{Bool}
        @test convert(SingleParam{String},str_single) isa SingleParam{String}
        @test convert(SingleParam{String},substr_single) isa SingleParam{String}
        @test convert(SingleParam{Float64},intbool) isa SingleParam{Float64}
        @test convert(SingleParam{Int},floatbool) isa SingleParam{Int}
        single_dual = Clapeyron.promote_model(Clapeyron.ForwardDiff.Dual{Nothing,Float64,2},floatbool)
        single_primal = Clapeyron.primalval(single_dual)
        @test single_primal.values == floatbool.values

        #SingleParam - indexing
        @test intbool["aa"] == 1
        intbool["bb"] = 2
        @test intbool["bb"] == 2
        @test_throws BoundsError intbool["cc"]
        @test_throws BoundsError setindex!(intbool,2,"cc")
        
        floatbool .= exp.(1.1 .+ floatbool)
        @test_throws InexactError convert(SingleParam{Int},floatbool)

        #pair param conversion and broadcasting
        floatbool = PairParam("float to bool",comps,[1.0 0.0;0.0 1.0])
        intbool = PairParam("int to bool",comps,[1 2;3 4])
        @test size(floatbool) == (2,2)
        @test convert(PairParam{Bool},floatbool) isa PairParam{Bool}
        @test convert(PairParam{Float64},intbool) isa PairParam{Float64}
        @test convert(PairParam{Int},floatbool) isa PairParam{Int}
        pair_dual = Clapeyron.promote_model(Clapeyron.ForwardDiff.Dual{Nothing,Float64,2},floatbool)
        pair_primal = Clapeyron.primalval(pair_dual)
        @test pair_primal.values == floatbool.values

        floatbool .= exp.(1.1 .+ floatbool)
        @test_throws InexactError convert(PairParam{Int},floatbool)
        #indexing for PairParam
        floatbool[1,2] = 1000
        @test floatbool[2,1] == 1000
        floatbool[1,2,false] = 4000
        @test floatbool[2,1] == 1000
        floatbool[1] = 1.2
        @test floatbool[1,1] == 1.2
        @test floatbool["aa","bb"] == 4000
        @test floatbool["aa"] == 1.2
        @test floatbool["aa","aa"] == 1.2
        @test floatbool["bb","aa"] == 1000
        @test_throws BoundsError floatbool["cc"]
        @test_throws BoundsError floatbool["cc","aa"]
        @test_throws BoundsError floatbool["aa","cc"]
        #pack vectors
        s1 = SingleParam("s1",comps,[1,2])
        s2 = SingleParam("s2",comps,[10,20])
        s = Clapeyron.pack_vectors(s1,s2)
        @test s[1] == [1,10]

        #assoc param
        ijab = [(1,1,1,2),(2,2,1,2),(2,2,1,3)]
        #assoc_sites = [["a1","a2"],["a1","a2"],["a1","a2","a3"]]
        
        float_vals = [1.,2.,3.]
        int_vals = [1,0,1]
        str_vals = ["v1","v2","v3"]

        assoc_vals_float = Clapeyron.Compressed4DMatrix(float_vals,ijab)
        assoc_vals_int = Clapeyron.Compressed4DMatrix(int_vals,ijab)
        assoc_vals_str = Clapeyron.Compressed4DMatrix(str_vals,ijab)
        assoc_vals_substr = Clapeyron.Compressed4DMatrix(chop.(str_vals),ijab)

        assoc_float = AssocParam("float",["comp1","comp2"],assoc_vals_float)
        assoc_int = AssocParam("int",["comp1","comp2"],assoc_vals_int)
        assoc_str = AssocParam("str",["comp1","comp2"],assoc_vals_str)
        assoc_substr = AssocParam("str",["comp1","comp2"],assoc_vals_substr)

        #AssocParam: indexing
        @test eltype(assoc_float) == Float64
        @test assoc_float[1] == [0.0 1.0; 1.0 0.0]
        @test assoc_float[2] == [0.0 2.0 3.0; 2.0 0.0 0.0; 3.0 0.0 0.0]
        @test iszero(prod(size(assoc_float[1,2])))

        #AssocParam: broadcasting
        assoc_float .= assoc_float .+ 2
        @test assoc_float.values.values == [3,4,5]
        assoc_float .= [1,1,1]
        @test assoc_float.values.values == [1,1,1]
        #AssocParam: Conversion
        @test convert(AssocParam{Bool},assoc_float) isa AssocParam{Bool}
        @test convert(AssocParam{Float64},assoc_int) isa AssocParam{Float64}
        @test convert(AssocParam{Int},assoc_int) isa AssocParam{Int}
        @test convert(AssocParam{String},assoc_substr) isa AssocParam{String}
        @test convert(AssocParam{String},assoc_str) isa AssocParam{String}
        assoc_dual = Clapeyron.promote_model(Clapeyron.ForwardDiff.Dual{Nothing,Float64,2},assoc_float)
        assoc_primal = Clapeyron.primalval(assoc_dual)
        @test assoc_primal.values.values == assoc_float.values.values

        assoc_float .= 1.2
        @test_throws InexactError convert(AssocParam{Int},assoc_float)
        @test_throws InexactError convert(AssocParam{Bool},assoc_float)
        assoc_float .= [3,4,5]
        #Compressed4DMatrix: conversion
        assoc_float_mat = [collect(assoc_float[i,j]) for (i,j) in Iterators.product(1:2,1:2)]
        assoc_float2 = AssocParam(assoc_float.name,assoc_float.components,assoc_float_mat)
        @test assoc_float.values.values == assoc_float2.values.values

        #single param: mixing rules with Real and Int type (#372)
        single_real = SingleParam("tes",["a","b","c"],Real[1,2,3])
        single_int = SingleParam("tes",["a","b","c"],[1,2,3])
        @test Clapeyron.sigma_LorentzBerthelot(single_real) isa PairParam
        @test Clapeyron.sigma_LorentzBerthelot(single_int) isa PairParam
        @test Clapeyron.lambda_LorentzBerthelot(single_real,0.5) isa PairParam
        @test Clapeyron.lambda_LorentzBerthelot(single_int,0.5) isa PairParam
    end

    
    @testset "GroupParam" begin
        # GC test, 3 comps, 4 groups
        components_gc = GroupParam(["test1", "test2", ("test3", ["grp1" => 2, "grp2" => 2, "grp3" => 3,"grp4" => 5])]; group_userlocations=filepath_gc)

        #Printing: GroupParam
        @test repr(components_gc) == "GroupParam[\"test1\" => [\"grp1\" => 1, \"grp2\" => 2], \"test2\" => [\"grp2\" => 1], \"test3\" => [\"grp1\" => 2, \"grp2\" => 2, \"grp3\" => 3, \"grp4\" => 5]]"
        @test repr("text/plain",components_gc) == "GroupParam(:test) with 3 components:\n \"test1\": \"grp1\" => 1, \"grp2\" => 2\n \"test2\": \"grp2\" => 1\n \"test3\": \"grp1\" => 2, \"grp2\" => 2, \"grp3\" => 3, \"grp4\" => 5"
        @test components_gc.grouptype == :test
        @test components_gc.components == ["test1", "test2", "test3"]
        @test components_gc.groups == [["grp1","grp2"],["grp2"],["grp1","grp2","grp3","grp4"]]
        @test components_gc.n_groups == [[1,2], [1], [2,2,3,5]]

        # Check that flattening of groups is correct.
        @test components_gc.flattenedgroups == ["grp1", "grp2", "grp3","grp4"]
        @test components_gc.n_flattenedgroups == [[1,2,0,0], [0,1,0,0], [2,2,3,5]]
        # Build param struct using the gc components above

        param_gc = getparams(components_gc; userlocations=filepath_param_gc)
        @test param_gc["param1"].values == [1, 2, 3, 4]
    end

    paramtable_file = ParamTable(:single,(species = ["sp1","sp2"],userparam = [2,10]))
    @testset "ParamTable" begin
        #reading external data, via ParamTable
        param_user = getparams(testspecies,userlocations = [paramtable_file],ignore_missing_singleparams=["userparam"])
        @test param_user["userparam"].values[1] === 2
    end

    @testset "direct csv parsing" begin
        #reading external data, via direct CSV parsing:
        csv_string = """Clapeyron Database File,
        in memory like parameters
        species,userparam,b,c
        sp1,1000,0.05,4
        sp2,,0.41,5
        """
        param_user2 = Clapeyron.getparams(["sp1","sp2"],userlocations = [csv_string],ignore_missing_singleparams=["userparam"])
        @test param_user2["userparam"].values[1] == 1000

        #@REPLACE keyword
        param_user3 = Clapeyron.getparams(["sp1","sp2"],userlocations = [paramtable_file, "@REPLACE/" * csv_string],ignore_missing_singleparams = ["userparam"])
        @test param_user3["userparam"].ismissingvalues[2] == true
    end

    @testset "named tuple parsing" begin

        model_nt = PCSAFT(["a1"],userlocations = (;
            Mw = [1.],
            epsilon = [2.],
            sigma = [3.],
            segment = [4.],
            k = fill(0.0,(1,1)),
            n_H = [1],
            n_e = [1],
            epsilon_assoc = Dict((("a1","e"),("a1","H")) => 1000.), 
            bondvol = Dict((("a1","e"),("a1","H")) => 0.001)))

        @test model_nt.params.Mw[1] == 1.
        @test model_nt.params.epsilon[1] == 2.
        @test model_nt.params.sigma[1] == 3e-10 #pcsaft multiplies sigma by 1e-10
        @test model_nt.params.segment[1] == 4.
        @test model_nt.params.epsilon_assoc.values.values[1] == 1000.
        @test model_nt.params.bondvol.values.values[1] == 0.001

        model_nt2 = PCSAFT(["a1"],userlocations = (;
            Mw = [1.],
            epsilon = [2.],
            sigma = [3.],
            segment = [4.],
            k = fill(0.0,(1,1)),
            epsilon_assoc = nothing, 
            bondvol = nothing))
        @test length(model_nt2.params.bondvol.values.values) == 0
    end

    @testset "misc database utils" begin
        #iszero
        @test Clapeyron._iszero(0.)
        @test Clapeyron._iszero(missing)
        @test Clapeyron._iszero("")

        #zero
        @test Clapeyron._zero(String) == ""
        #@test Clapeyron._zero(Missing) == "" #need to check if we actually use this
        @test Clapeyron._zero(Int) == 0

        #defaultmissing
        @test Clapeyron.defaultmissing(["a","b","",missing])[2] == Bool[0, 0, 0, 1]
        @test Clapeyron.defaultmissing(["a","b","","d"])[2] == Bool[0, 0, 0, 0]
        @test Clapeyron.defaultmissing([true,false,missing])[2] == Bool[0, 0, 1]
        @test all(Clapeyron.defaultmissing(fill(missing,4))[2])
        @test Clapeyron.defaultmissing(["a",1,missing])[1] == ["a","1",""]
        @test Clapeyron.defaultmissing(["a",1,missing])[2] == Bool[0, 0, 1]

        c1 = Clapeyron.cas("water")
        @test c1[1] == "7732-18-5"
    end
end