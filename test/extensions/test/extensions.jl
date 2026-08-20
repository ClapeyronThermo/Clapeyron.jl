Clapeyron.use_superancillaries!(false)

@testset "Superancillaries.jl" begin
        
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

@testset "MultiComponentFlash.jl Algorithm" begin

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