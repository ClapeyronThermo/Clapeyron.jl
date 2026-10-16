const  eos_ = [cPR] # list of EoS models to test

@testset verbose = true "Carnot.jl" begin

fluids_test = ["dimethyl ether","isobutane","Butene","Neopentane","isopentane",
            "pentane","hexane","acetone","cyclopentane","heptane","cyclohexane","benzene","propane"]

    solverparams_ad = ThermoCycleParameters(
    autodiff = true,
    xtol = 1e-8,
    ftol = 1e-8,
    restart_TOL = 1e-4,
    max_iters = 1000,
    N = 30,
    x0_init = :default
)

solverparams_fd = ThermoCycleParameters(
    autodiff = false,
    xtol = 1e-8,
    ftol = 1e-8,
    restart_TOL = 1e-4,
    max_iters = 1000,
    N = 30,
    x0_init = :default
)

@testset "Isentropic Compression - Single Component" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z1 = [1.0]; p1 = 101325;
        T1 = dew_temperature(fluid_model,p1,z1)[1] + 50;  p2 = 101325*2;
        h1 = enthalpy(fluid_model,p1,T1,z1);
        s1 = Clapeyron.entropy(fluid_model,p1,T1,z1)
        h2= Carnot.isentropic_compressor(p1,p2,1.0,h1,z1,fluid_model)
        s2 = Clapeyron.PH.entropy(fluid_model,p2,h2,z1)
        @test isapprox(s1,s2,atol=1e-5)
    end
end

@testset "Isentropic pump - Single Component" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z1 = [1.0]; p1 = 101325;
        T1 = bubble_temperature(fluid_model,p1,z1)[1] - 20;  p2 = 101325*2;
        h1 = enthalpy(fluid_model,p1,T1,z1);
        s1 = Clapeyron.entropy(fluid_model,p1,T1,z1)
        h2= Carnot.isentropic_pump(p1,p2,1.0,h1,z1,fluid_model)
        s2 = Clapeyron.PH.entropy(fluid_model,p2,h2,z1)
        @test isapprox(s1,s2,atol=1e-5)
    end
end

@testset "Isentropic Expansion - Single Component" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z1 = [1.0]
        p1 = 101325*5;T1 = saturation_temperature(fluid_model,p1)[1] + 100.0; p2 = 101325;
        h1 = enthalpy(fluid_model,p1,T1,z1);
        s1 = Clapeyron.entropy(fluid_model,p1,T1,z1)
        h2= Carnot.isentropic_expander(p1,p2,1.0,h1,z1,fluid_model)
        s2 = Clapeyron.PH.entropy(fluid_model,p2,h2,z1)
        @test isapprox(s1,s2,atol=1e-5)
    end
end

@testset "ORC - mixture AD" begin
    fluid = cPR(["propane","butane"],idealmodel = ReidIdeal);
    T_evap_out = 310
    z = [0.6157894736842106, 9.384210526315789]
    ΔT_sh = 11.473684210526315 
    _orc_ = ORC(fluid=fluid, z=z, T_evap_in=330, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
    sol = solve(_orc_,solverparams_ad)
    @test Carnot.norm(sol.residuals) < 1e-3
    @test sol.soltype == :subcritical
end

@testset "ORC - pure AD" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z = [1.0]
        T_evap_out = 310
        ΔT_sh = 5.0 
        _orc_ = ORC(fluid=fluid_model, z=z, T_evap_in=330, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
        sol = solve(_orc_,solverparams_ad)
        @test Carnot.norm(sol.residuals) < 1e-3
        @test sol.soltype == :subcritical
    end
end

@testset "ORC - pure AD hard" begin
    for fluid in fluids_test
        fluid_model = cPR(fluid,idealmodel = ReidIdeal)
        z = [1.0]
        T_evap_out = 310
        ΔT_sh = 5.0 
        _orc_ = ORC(fluid=fluid_model, z=z, T_evap_in=330, T_evap_out=T_evap_out, T_cond_in=270, T_cond_out=280, η_pump=0.75, η_expander=0.75, pp_evap=3, pp_cond=3, ΔT_sc=2.0, ΔT_sh=ΔT_sh)
        sol = solve(_orc_,solverparams_ad)
        @test Carnot.norm(sol.residuals) < 1e-3
        @test sol.soltype == :subcritical
    end
end

@testset "issue : #20" begin
    fluid = cPR(["acetone", "isopentane"],idealmodel= ReidIdeal); z = [0.1, 9.9];
    orc = ORC(
        fluid = fluid, z= z,
        T_evap_in = 383.15, T_evap_out = 363.15, 
        T_cond_in = 290.0, T_cond_out = 310.0,
        η_pump = 0.75, η_expander = 0.75,
        pp_evap = 5.0, pp_cond = 5.0,
        ΔT_sc = 16.174645986149887, ΔT_sh = 6.851340038880649
    )
    solver_params = ThermoCycleParameters()
    sol = solve(orc,solver_params)
    @test Carnot.norm(sol.residuals) < 1e-3
    @test sol.soltype == :subcritical 
    @test sol.x[1]>=sol.x[2]
end

end