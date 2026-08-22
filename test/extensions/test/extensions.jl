Clapeyron.use_superancillaries!(false)

@testset "Superancillaries.jl Extension" begin
        
    pc = PCSAFT("eicosane")
    pc2 = pharmaPCSAFT("oxygen")
    cubic = tcPR(["water"])
    
    crit_pc = crit_pure(pc)
    sat_cubic = saturation_pressure(cubic,373.15)
    sat_pc2 = saturation_pressure(pc2,150.0)
    
    Clapeyron.use_superancillaries!(true)
    
    crit_sa_pc = crit_pure(pc)
    sat_sa_cubic = saturation_pressure(cubic,373.15)
    sat_sa_pc2 = saturation_pressure(pc2,150.0)

    @test crit_pc[1] ≈ crit_sa_pc[1] rtol = 1e-6
    @test crit_pc[3] ≈ crit_sa_pc[3] rtol = 1e-6
    
    @test sat_cubic[1] ≈ sat_sa_cubic[1] rtol = 1e-6
    @test sat_cubic[2] ≈ sat_sa_cubic[2] rtol = 1e-6
    @test sat_cubic[3] ≈ sat_sa_cubic[3] rtol = 1e-6
    
    @test sat_pc2[1] ≈ sat_sa_pc2[1] rtol = 1e-6
    @test sat_pc2[2] ≈ sat_sa_pc2[2] rtol = 1e-6
    @test sat_pc2[3] ≈ sat_sa_pc2[3] rtol = 1e-6
end

@testset "MultiComponentFlash.jl Extension" begin
    system = PCSAFT(["water","cyclohexane","propane"])
    T = 298.15
    p = 1e5
    z = [0.333, 0.333, 0.334]
    #two-phase test, using Clapeyron api
    mcf = MCFlashJL()
    @test Clapeyron.numphases(Clapeyron.tp_flash2(system, p, T, z, mcf)) == 2
    #vapour test, using MCF api
    cond = (p = 5e6, T = 303.15, z = [0.4, 0.6])
    vapour_model = PR78(["hydrogen", "methane"])
    vapour_res = MultiComponentFlash.flashed_mixture_2ph(vapour_model,cond)
    @test vapour_res.state == MultiComponentFlash.single_phase_v
    @test vapour_res.vapor.Z ≈ 0.9672507136048648 rtol = 1e-6

    #liquid test,using MCF api
    liquid_model = cPR(["octane","nonane"])
    cond = (p = 5e7, T = 303.15, z = [0.4, 0.6])
    liquid_res = MultiComponentFlash.flashed_mixture_2ph(liquid_model,cond)
    @test liquid_res.state == MultiComponentFlash.single_phase_l
    @test liquid_res.liquid.Z ≈ 3.458550315299117 rtol = 1e-6
end


tdb(calc::Glenn.Calculator) = calc.db
tdb(db::Glenn.ThermoDB) = db

@testset "Glenn.jl Extension" begin
    calc = Glenn.Calculator()

    @testset "Model creation" begin
        model_O2 = Clapeyron.GlennJL(calc, "O2"; R0=8.31446261815324)
        @test model_O2 isa GlennJL
        @test length(model_O2) == 1
        @test model_O2.components[1] == "O2"
        @test model_O2.species_info[1].phase == "gas"

        model_mix = Clapeyron.GlennJL(calc, ["N2", "O2"])
        @test length(model_mix) == 2
        @test model_mix.components == ["N2", "O2"]
        @test all(phase -> phase == "gas", [info.phase for info in model_mix.species_info])

        spN2 = only(Glenn.find_species(tdb(calc), "N2", exact_match=true))
        spO2 = only(Glenn.find_species(tdb(calc), "O2", exact_match=true))
        model_idx = Clapeyron.GlennJL(calc, [spN2.id, spO2.id])
        @test model_idx.components == ["N2", "O2"]

        model_info = Clapeyron.GlennJL(calc, [spN2, spO2])
        @test model_info.components == ["N2", "O2"]

        @test_throws ErrorException Clapeyron.GlennJL(calc, [1.0, 2.0])
    end

    @testset "Molecular weight" begin
        model_O2 = Clapeyron.GlennJL(calc, "O2")
        @test Clapeyron.molecular_weight(model_O2, [1.0]) ≈ 0.0319988
        @test Clapeyron.mw(model_O2)[1] ≈ 31.9988

        model_mix = Clapeyron.GlennJL(calc, ["N2", "O2"])
        @test Clapeyron.molecular_weight(model_mix, [0.5, 0.5]) ≈ 0.5*28.0134e-3 + 0.5*31.9988e-3
        @test Clapeyron.mw(model_mix) == [28.0134, 31.9988]
    end

    @testset "Pure species properties" begin
        model_O2 = Clapeyron.GlennJL(calc, "O2")
        T = 300.0
        sp_info = model_O2.species_info[1]
        sp_data = Glenn.get_species_data(calc.db, sp_info.id)
        interval = sp_data["intervals"][1]
        coeff = interval.coefficients

        cp_glenn = Glenn.calculate_cp(model_O2, T)
        cp_direct = Glenn.ThermoDatabase.calculate_cp(coeff, T) * model_O2.R0[1]
        @test cp_glenn ≈ cp_direct atol=1e-9

        h_glenn = Glenn.calculate_h(model_O2, T)
        s_glenn = Glenn.calculate_s(model_O2, T)
        prop = Glenn.calculate_properties(model_O2, T)
        @test prop.temperature == T
        @test prop.cp ≈ cp_glenn
        @test prop.h_relative ≈ h_glenn * T * Rgas(model_O2)
        @test prop.s ≈ s_glenn * Rgas(model_O2)
        #@test prop.temp_min == interval.temp_min
        #@test prop.temp_max == interval.temp_max
        @test prop.species_name == "O2"
        @test prop.phase == "gas"
    end

    @testset "Mixture properties" begin
        model_mix = Clapeyron.GlennJL(calc, ["N2", "O2"])
        T = 500.0
        z = [0.8, 0.2]

        cp_mix = Glenn.calculate_cp(model_mix, T, z)
        cp_N2 = Glenn.calculate_cp(Clapeyron.GlennJL(calc, "N2"), T)
        cp_O2 = Glenn.calculate_cp(Clapeyron.GlennJL(calc, "O2"), T)
        cp_expected = (z[1]*cp_N2 + z[2]*cp_O2) / sum(z)
        @test cp_mix ≈ cp_expected atol=1e-9

        h_mix = Glenn.calculate_h(model_mix, T, z)
        h_N2 = Glenn.calculate_h(Clapeyron.GlennJL(calc, "N2"), T)
        h_O2 = Glenn.calculate_h(Clapeyron.GlennJL(calc, "O2"), T)
        h_expected = (z[1]*h_N2 + z[2]*h_O2) / sum(z)
        @test h_mix ≈ h_expected atol=1e-9

        s_mix = Glenn.calculate_s(model_mix, T, z)
        s_N2 = Glenn.calculate_s(Clapeyron.GlennJL(calc, "N2"), T)
        s_O2 = Glenn.calculate_s(Clapeyron.GlennJL(calc, "O2"), T)
        s_expected = (z[1]*s_N2 + z[2]*s_O2) / sum(z)
        @test s_mix ≈ s_expected atol=1e-9

        hf_mix = Glenn.calculate_formation_enthalpy(model_mix, z)
        hf_N2 = model_mix.species_info[1].heat_of_formation_298K
        hf_O2 = model_mix.species_info[2].heat_of_formation_298K
        hf_expected = (z[1]*hf_N2 + z[2]*hf_O2) / sum(z)
        @test hf_mix ≈ hf_expected

        T1, T2 = 300.0, 1000.0
        Δh = Glenn.calculate_enthalpy_change(model_mix, T1, T2, z)
        h1 = Glenn.calculate_h(model_mix, T1, z) * T1 * Rgas(model_mix)
        h2 = Glenn.calculate_h(model_mix, T2, z) * T2 * Rgas(model_mix)
        @test Δh ≈ h2 - h1

        prop_mix = Glenn.calculate_properties(model_mix, T, z)
        @test prop_mix.species_name == "mixture"
        @test prop_mix.phase == "gas"
        @test prop_mix.cp ≈ cp_mix
        @test prop_mix.h_relative ≈ h_mix * T * Rgas(model_mix)
        @test prop_mix.s ≈ s_mix * Rgas(model_mix)
        @test isnan(prop_mix.temp_min)
        @test isnan(prop_mix.temp_max)
    end

    @testset "Clapeyron ideal gas functions" begin
        model_O2 = Clapeyron.GlennJL(calc, "O2")
        V = 1.0
        T = 500.0
        z = [1.0]
        @test Clapeyron.a_ideal(model_O2, V, T, z) ≈ -29.260376033957428 rtol = 1e-6

        d2f = Clapeyron.∂²f∂T²(model_O2, V, T, z)
        @test d2f < 0
        model_mix = Clapeyron.GlennJL(calc, ["N2", "O2"])
        z_mix = [0.7, 0.3]
        @test Clapeyron.a_ideal(model_mix, V, T, z_mix) ≈ -28.725485094123094 rtol = 1e-6
        d2f_mix = Clapeyron.∂²f∂T²(model_mix, V, T, z_mix)
        @test d2f_mix < 0
    end

    @testset "get_properties_range" begin
        model_O2 = Clapeyron.GlennJL(calc, "O2")
        T_range = [300.0, 400.0, 500.0]
        props = Glenn.get_properties_range(model_O2, T_range)
        @test length(props) == 3
        for prop in props
            @test prop isa ThermoProperties
        end
    end
end