ENV["JULIA_TEST_FAILFAST"] = "false"
using Test
t1 = @elapsed using Clapeyron
using CoolProp #CoolProp ext
using Unitful #Unitful ext
using MultiComponentFlash: MultiComponentFlash
using Clapeyron.LinearAlgebra
using Clapeyron.StaticArrays
using Clapeyron: has_sites,has_groups

#=

Modify this constant to true to run all tests in all workers

=#

ALL_TESTS = false
include("utils.jl")
@info "Loading Clapeyron took $(round(t1,digits = 2)) seconds"
@info "Coolprop: $(Clapeyron.is_coolprop_loaded())"
#Disable showing citations
ENV["CLAPEYRON_SHOW_REFERENCES"] = "FALSE"

#fix to current tests
function GERG2008(components;verbose = false,reference_state = nothing)
    return MultiFluid(components;
    mixing = AsymmetricMixing,
    departure = EmpiricDeparture,
    pure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/pures"],
    mixing_userlocations  = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/mixing/GERG2008_mixing_unlike.csv"],
    departure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/departure/GERG2008_departure_unlike.csv"],
    reference_state = reference_state,
    coolprop_userlocations = false,
    verbose = verbose,
    Rgas = Clapeyron.R̄)
end

function test_gibbs_duhem(model,V,T,z;rtol = 1e-14)
    for i in (2.0,3.0,5.0,7.0,11.0)
        a_res₀ = Clapeyron.a_res(model,V,T,z)
        @test a_res₀ ≈ Clapeyron.a_res(model,i*V,T,i*z) rtol = rtol
    end
    pures = split_model(model)
    x_pure = zeros(length(model))
    for n in 1:length(model) 
        for i in (2.0,3.0,5.0,7.0,11.0)
            x_pure[n] = i
            @test Clapeyron.a_res(model,i*V,T,x_pure) ≈ Clapeyron.a_res(pures[n],i*V,T,Clapeyron.SVector(i)) rtol = rtol
            x_pure .= 0
        end
    end
end

function test_volume(model,p,T,z = Clapeyron.SA[1.0],rtol = 1e-8)
    v = volume(model,p,T,z)
    @test p ≈ Clapeyron.pressure(model,v,T,z) rtol = rtol
end

function test_scales(model,z)
    z3 = 3 .* z
    z7 = 7 .* z
    T0 = 300.0
    @test lb_volume(model,T0,z3) ≈ 3*lb_volume(model,z)
    @test lb_volume(model,T0,z7) ≈ 7*lb_volume(model,z)
    @test T_scale(model,z3) ≈ T_scale(model,z)
    @test p_scale(model,z3) ≈ p_scale(model,z)
end
#=
include_distributed distributes the test load among all workers
=#
display(Test.detect_ambiguities(Clapeyron))
include_distributed("test_database.jl",4)
include_distributed("test_solvers.jl",4)
include_distributed("test_differentials.jl",4)
include_distributed("test_misc.jl",4)
include_distributed("test_models_saft_pc.jl",1)
include_distributed("test_models_cubic.jl",3)
include_distributed("test_models_saft_others.jl",3)
include_distributed("test_models_others.jl",2)
include_distributed("test_models_saft_vr.jl",1)
include_distributed("test_models_electrolytes.jl",1)
include_distributed("test_methods_eos.jl",4)
include_distributed("test_methods_api.jl",2)
include_distributed("test_methods_api_flash.jl",3)
include_distributed("test_methods_electrolytes.jl",1)
include_distributed("test_estimation.jl",1)
include_distributed("test_issues.jl",1)

