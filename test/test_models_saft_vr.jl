GC.gc()

@testset "SAFT-VR-Mie Models" begin
    @printline
    let T = 298.15, V = 1e-4,z1 = Clapeyron.SA[1.0],z = [0.5,0.5],z3 = [0.333, 0.333,0.333];
    @testset "SAFTVRMie" begin
        system = SAFTVRMie(["methanol", "water"])
        @test Clapeyron.a_mono(system, V, T, z) ≈ -0.9729139704318698 rtol = 1e-6
        _a_chain = Clapeyron.a_chain(system, V, T, z)
        _a_disp  = Clapeyron.a_disp(system, V, T, z)
        @test _a_chain ≈ -0.028347378889242814 rtol = 1e-6
        @test Clapeyron.a_dispchain(system,V,T,z) - _a_chain ≈ _a_disp rtol = 1e-6
        @test Clapeyron.a_assoc(system, V, T, z) ≈ -4.18080707238976 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        test_recombine(system)
        GC.gc()
    end

    @testset "SAFTVRMieGV" begin
        system = SAFTVRMieGV(["benzene","acetone"])
        V_GV = 8e-5
        test_gibbs_duhem(system,V_GV,T,z)
        test_recombine(system)
        
        
        @test Clapeyron.a_mp(system, V_GV, T, z) ≈ -0.7521858819355216 rtol = 1e-6
        
        system2 = SAFTVRMieGV(["carbon dioxide","benzene"])
        @test Clapeyron.a_mp(system2, V_GV, T, z) ≈ -0.33476290652200424 rtol = 1e-6

        system3 = SAFTVRMieGV(["acetone","diethyl ether","ethyl acetate"])
        @test Clapeyron.a_mp(system3, V_GV, T, [0.3,0.3,0.4]) ≈ -0.8088095725122674 rtol = 1e-6
        GC.gc()
    end

    @testset "SAFTVRQMie" begin
        system = SAFTVRQMie(["helium"])
        @test Clapeyron.a_mono(system, V, T, z1) ≈ 0.12286776703976324 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z1)
        GC.gc()
    end

    @testset "SAFTVRSMie" begin
        system = SAFTVRSMie(["carbon dioxide"])
        V_sol = 3e-5
        @test Clapeyron.a_mono(system, V_sol, T, z1) ≈ 0.43643302846919896 rtol = 1e-6
        @test Clapeyron.a_chain(system, V_sol, T, z1) ≈ -0.4261294644079463 rtol = 1e-6
        test_gibbs_duhem(system,V_sol,T,z1,rtol = 1e-12)
        GC.gc()
    end

    @testset "SAFTVRMie15" begin
        v15 = 1/(1000*1000*1.011/18.015)
        T15 = 290.0
        vr15 = SAFTVRMie15("water")
        #Dufal, table 4, 290K, f_OH(free) = 0.089
        @test Clapeyron.X(vr15,v15,T15,z1)[1][1] ≈ 0.08922902098124778 rtol = 1e-6
    end

    @testset "SAFTgammaMie" begin
        system = SAFTgammaMie(["methanol","butane"])
        V_γMie = exp10(-3.5)
        @test Clapeyron.a_mono(system, V_γMie, T, z) ≈ -1.0400249396482548 rtol = 1e-6
        @test Clapeyron.a_chain(system, V_γMie, T, z) ≈ -0.07550931466871749 rtol = 1e-6
        @test Clapeyron.a_assoc(system, V_γMie, T, z) ≈ -0.8205840455850311 rtol = 1e-6
        test_gibbs_duhem(system,V,T,z)
        test_scales(system)
        test_recombine(system)
        test_repr(system,str = ["\"methanol\": \"CH3OH\" => 1"])
        GC.gc()

        #=
        #171
        This is a problem that occurs in an intersection between split_model and cross-association sites
        a single component model created from scratch don't have any cross association sites,
        but a single component model created from split_model does have those sites.
        we neet to check if the results of both are equal.
        =#
        m171 = SAFTgammaMie(["water","acetone"])
        m171_split = split_model(m171)[2]
        m171_pure = SAFTgammaMie(["acetone"])
        res_pure = Clapeyron.eos(m171_pure,1.013e6,298.15) #works
        res_split = Clapeyron.eos(m171_split,1.013e6,298.15) #should work
        @test res_pure ≈ res_split
    end

    @testset "SAFTgammaMie - DM error 1" begin
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

    @testset "SAFTgammaMie - #145" begin
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

    @testset "SAFTgammaMie - #262" begin
        #=
        SAFT-gamma Mie
        Ensure ij,ab == ji,ba, and that the group-comp is correct for assoc
        =#
        model = SAFTgammaMie(["water","ethanol"])
        @test model.params.epsilon_assoc.values[1,2][1,2] == model.params.epsilon_assoc.values[2,1][2,1]
        @test model.params.epsilon_assoc.values[1,2][1,2] == model.vrmodel.params.epsilon_assoc.values[1,2][1,2]
    end

    @testset "SAFTgammaMie_custom_name_with_MonomerIdeal" begin
        species = [("test_propane", ["CH3"=>2, "CH2"=>1])]
        system = SAFTgammaMie(species; idealmodel = MonomerIdeal)
        @test !isnothing(system)
        GC.gc()

        #GC model, GC ideal model
        system2 = SAFTgammaMie(species; idealmodel = JobackIdeal)
        @test !isnothing(system)
    end
    
    @testset "structSAFTgammaMie" begin
        species = [("ethanol",["CH3"=>1,"CH2OH"=>1],[("CH3","CH2OH")=>1]),
                   ("octane",["CH3"=>2,"CH2"=>6],[("CH3","CH2")=>2,("CH2","CH2")=>5])]

        system = structSAFTgammaMie(species)
        V_γMie = exp10(-3.5)
        @test Clapeyron.a_chain(system, V_γMie, T, z) ≈ -0.11160851237651681 rtol = 1e-6
        test_gibbs_duhem(system,V_γMie,T,z,rtol = 1e-12)
        GC.gc()
    end
    end
    @printline
end
