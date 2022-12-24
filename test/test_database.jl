using Clapeyron, Test, LinearAlgebra

@testset "database_lookup" begin

    @test Clapeyron.normalisestring("COO-") == "coo-" #before, it was coo, making a collision on Electrolyte SAFTgammaMie

    params1 = Clapeyron.getparams(["water", "methanol"], ["SAFT/PCSAFT"],return_sites=false)
    @test haskey(params1, "sigma")

    @printline
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
    @test_throws ErrorException Clapeyron.compile_params(testspecies,allparams,allnotfoundparams,opts) #generate ClapeyronParams
    result, allcomponentsites = Clapeyron.compile_params(testspecies,allparams,allnotfoundparams,opts2) #generate ClapeyronParams

    @test allcomponentsites == [[],
                                                                     [],
                                                                     ["e", "e2", "H"],
                                                                     ["e", "H"],
                                                                     ["e", "e2", "H"]]

    # Check that it throws an error if ignore_missing_singleparams is not set to true.
    @test_throws ErrorException Clapeyron.getparams(testspecies; userlocations=filepath_normal,return_sites = false)

    params = Clapeyron.getparams(testspecies; userlocations=filepath_normal, ignore_missing_singleparams=["emptyparam","missingparam"],return_sites = false)
    sites = Clapeyron.SiteParam(params["intparam"].components)
    #Printing: SingleParam
    @test repr(params["intparam"]) == "SingleParam{Int64}(\"intparam\")[\"sp1\", \"sp2\", \"sp3\", \"sp4\", \"sp5\"]"
    @test repr("text/plain",params["intparam"]) == "SingleParam{Int64}(\"intparam\") with 5 components:\n \"sp1\" => 6\n \"sp2\" => 2\n \"sp3\" => 7\n \"sp4\" => 4\n \"sp5\" => 5"
    #Printing: PairParam
    @test repr(params["overwriteparam"]) == "PairParam{Float64}(\"overwriteparam\")[\"sp1\", \"sp2\", \"sp3\", \"sp4\", \"sp5\"]"
    @test repr("text/plain",params["overwriteparam"]) == "5×5 PairParam{Float64}([\"sp1\", \"sp2\", \"sp3\", \"sp4\", \"sp5\"]) with values:\n 1.6  4.0  0.0  0.0  0.0\n 4.0  1.2  0.0  3.0  0.0\n 0.0  0.0  1.3  2.0  0.0\n 0.0  3.0  2.0  1.4  0.0\n 0.0  0.0  0.0  0.0  1.5"
    #Printing: AssocParam
    @test repr(params["overwriteassocparam"]) == "AssocParam{String}(\"overwriteassocparam\")[\"val1\", \"val8\", \"val5\", \"val4\", \"val7\", \"val6\", \"val3\", \"42\"]"

    @test repr("text/plain",params["overwriteassocparam"]) == "AssocParam{String}[\"sp1\", \"sp2\", \"sp3\", \"sp4\", \"sp5\"]) with 8 values:\n(\"sp3\", \"e\") >=< (\"sp3\", \"H\"): val1\n(\"sp3\", \"e2\") >=< (\"sp3\", \"H\"): val8\n(\"sp3\", \"H\") >=< (\"sp4\", \"e\"): val5\n(\"sp3\", \"H\") >=< (\"sp4\", \"H\"): val4\n(\"sp4\", \"e\") >=< (\"sp5\", \"H\"): val7\n(\"sp4\", \"H\") >=< (\"sp5\", \"e2\"): val6\n(\"sp5\", \"e\") >=< (\"sp5\", \"e\"): val3\n(\"sp5\", \"e\") >=< (\"sp5\", \"H\"): 42"
    #Printing: SiteParam
    @test repr(sites) == "SiteParam[\"sp1\" => [], \"sp2\" => [], \"sp3\" => [], \"sp4\" => [], \"sp5\" => []]"
    @test repr("text/plain",sites) == "SiteParam with 5 components:\n \"sp1\": (no sites)\n \"sp2\": (no sites)\n \"sp3\": (no sites)\n \"sp4\": (no sites)\n \"sp5\": (no sites)"
    # Check that all the types are correct.
    @test typeof(params["intparam"]) <: Clapeyron.SingleParam{Int}
    @test typeof(params["doubleparam"]) <: Clapeyron.SingleParam{Float64}

    #@test typeof(params["boolparam"]) <: Clapeyron.SingleParam{Bool}
    #we can now convert directly, no need to parse diferently

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

    assoc_param_values =
    [[Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,3)          ]  [Array{Int64}(undef,0,2)]  [Array{Int64}(undef,0,3)       ]
     [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,3)          ]  [Array{Int64}(undef,0,2)]  [Array{Int64}(undef,0,3)       ]
     [Array{Int64}(undef,3,0)]  [Array{Int64}(undef,3,0)]  [[0 0 2000; 0 0 2700; 2000 2700 0]]  [[0 0; 0 0; 2400 2300]  ]  [[0 0 0; 0 0 0; 0 0 0]         ]
     [Array{Int64}(undef,2,0)]  [Array{Int64}(undef,2,0)]  [[0 0 2400; 0 0 2300]             ]  [[0 0; 0 0]             ]  [[0 0 2600; 0 2500 0]          ]
     [Array{Int64}(undef,3,0)]  [Array{Int64}(undef,3,0)]  [[0 0 0; 0 0 0; 0 0 0]            ]  [[0 0; 0 2500; 2600 0]  ]  [[2200 0 2100; 0 0 0; 2100 0 0]]]

    for i ∈ 1:5
        for j ∈ 1:5
            valij = assoc_param_values[i,j]
            valji = assoc_param_values[j,i]
            s1,s2 = size(valij)
            testij = params["assocparam"].values[i,j]
            if iszero(s1*s2)
                continue # a compressed assoc matrix can hold more information that is discarded normally
            else
                for a ∈ 1:s1
                    for b ∈ 1:s2
                        #this is to account for the symmetric nature of the CompressedAssoc4DMatrix
                        #where as the original input matrix is asymmetric
                        val_ij_ab = valij[a,b]
                        val_ji_ba = valji[b,a]
                        if !iszero(val_ij_ab)
                            val = val_ij_ab
                        else
                            val = val_ji_ba
                        end
                        @test testij[a,b] == val
                    end
                end
            end
        end
    end
  #  @test params["assocparam"].values[i,j] .== assoc_param_values[i,j] for (i,j) in zip()

    # Clashing headers between association and non-association parameters are not allowed

    @test_throws ErrorException Clapeyron.getparams(testspecies; userlocations=filepath_clashingheaders)

    # If parameter is not tagged as ignore_missing_singleparams, incomplete diagonal will throw an error
    #ignore_missing_singleparams=["asymmetricpair"]
    @test_throws ErrorException Clapeyron.getparams(testspecies; userlocations=filepath_asymmetry,return_sites = false)

    asymmetricparams = Clapeyron.getparams(testspecies; userlocations=filepath_asymmetry, asymmetricparams=["asymmetricpair", "asymmetricassoc"], ignore_missing_singleparams=["asymmetricpair"],return_sites = false)

    asymmetricparams["asymmetricpair"].values == [0.06  0.04  0.0  0.0   0.0
                                                  0.05  0.0   0.0  0.0   0.0
                                                  0.0   0.0   0.0  0.02  0.0
                                                  0.0   0.03  0.0  0.0   0.0
                                                  0.0   0.0   0.0  0.0   0.0]
    asymmetricparams["asymmetricassoc"].values ==
    [[Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,3)       ]  [Array{Int64}(undef,0,2) ]  [Array{Int64}(undef,0,3)    ]
     [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,0)]  [Array{Int64}(undef,0,3)       ]  [Array{Int64}(undef,0,2) ]  [Array{Int64}(undef,0,3)    ]
     [Array{Int64}(undef,3,0)]  [Array{Int64}(undef,3,0)]  [[0 0 2800; 0 0 0; 2000 2700 0]]  [[0 2900; 0 0; 2400 2300]]  [[0 0 0; 0 0 0; 0 0 0]      ]
     [Array{Int64}(undef,2,0)]  [Array{Int64}(undef,2,0)]  [[0 0 0; 0 0 0]                ]  [[0 0; 0 0]              ]  [[0 0 0; 0 2500 0]          ]
     [Array{Int64}(undef,3,0)]  [Array{Int64}(undef,3,0)]  [[0 0 0; 0 0 0; 0 0 0]         ]  [[0 0; 0 0; 2600 0]      ]  [[2200 0 0; 0 0 0; 2100 0 0]]]

    # Also, since a non-missing value exists on the diagonal of "asymmetricpair",
    # and the diagonal contains missing values, it should throw an error
    # when parameter is not in ignore_missing_singleparams
    @test_throws ErrorException Clapeyron.getparams(testspecies; userlocations=filepath_asymmetry, asymmetricparams=["asymmetricpair", "asymmetricassoc"])

    # Testing for multiple identifiers
    multiple_identifiers = Clapeyron.getparams(testspecies; userlocations=filepath_multiple_identifiers, return_sites=false)
    @test multiple_identifiers["param_single"].values == [100, 200, 200, 300, 300]
    @test multiple_identifiers["param_pair"].values   == [1000  4000     0     0     0
                                                         4000  2000     0  6000     0
                                                            0     0  2000  5000     0
                                                            0  6000  5000  3000     0
                                                            0     0     0     0  3000]
    @test multiple_identifiers["param_assoc"].values[3,2][1,1] == 10000 &&
            multiple_identifiers["param_assoc"].values[5,1][1,1] == 20000


    #single param conversion and broadcasting
    comps = ["aa","bb"]
    floatbool = SingleParam("float to bool",comps,[1.0,0.0])
    intbool = SingleParam("int to bool",comps,[1,0])
    @test size(floatbool) == (2,)
    @test convert(SingleParam{Bool},floatbool) isa SingleParam{Bool}
    @test convert(SingleParam{Float64},intbool) isa SingleParam{Float64}
    @test convert(SingleParam{Int},floatbool) isa SingleParam{Int}
    floatbool .= exp.(1.1 .+ floatbool)
    @test_throws AssertionError convert(SingleParam{Int},floatbool)

    #pair param conversion and broadcasting
    floatbool = PairParam("float to bool",comps,[1.0 0.0;0.0 1.0])
    intbool = PairParam("int to bool",comps,[1 2;3 4])
    @test size(floatbool) == (2,2)
    @test convert(PairParam{Bool},floatbool) isa PairParam{Bool}
    @test convert(PairParam{Float64},intbool) isa PairParam{Float64}
    @test convert(PairParam{Int},floatbool) isa PairParam{Int}
    floatbool .= exp.(1.1 .+ floatbool)
    @test_throws AssertionError convert(PairParam{Int},floatbool)
    #indexing for PairParam
    floatbool[1,2] = 1000
    @test floatbool[2,1] == 1000
    floatbool[1,2,false] = 4000
    @test floatbool[2,1] == 1000
    floatbool[1] = 1.2
    @test floatbool[1,1] == 1.2

    #pack vectors
    s1 = SingleParam("s1",comps,[1,2])
    s2 = SingleParam("s2",comps,[10,20])
    s = Clapeyron.pack_vectors(s1,s2)
    @test s[1] == [1,10]

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
    @test components_gc.n_flattenedgroups == [[1,2,0,0], [0,1,0,0 ], [2,2,3,5]]
    # Build param struct using the gc components above

    param_gc = getparams(components_gc; userlocations=filepath_param_gc)
    @test param_gc["param1"].values == [1, 2, 3, 4]

    #reading external data, via ParamTable
    file = ParamTable(:single,(species = ["sp1","sp2"],userparam = [2,10]))
    param_user = getparams(testspecies,userlocations = [file],ignore_missing_singleparams=["userparam"])
    @test param_user["userparam"].values[1] === 2

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
    param_user3 = Clapeyron.getparams(["sp1","sp2"],userlocations = [file, "@REPLACE/" * csv_string],ignore_missing_singleparams = ["userparam"])
    @test param_user3["userparam"].ismissingvalues[2] == true
end

