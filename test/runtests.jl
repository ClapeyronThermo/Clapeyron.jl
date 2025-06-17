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

function test_volume(model,p,T,z = Clapeyron.SA[1.0];rtol = 1e-8,phase = :unknown)
    v = volume(model,p,T,z)
    @test p ≈ Clapeyron.pressure(model,v,T,z) rtol = rtol
end

function test_scales(model,T0 = 300.0)
    n = length(model)
    z = ones(n) ./ n
    z3 = 3 .* z
    z7 = 7 .* z
    @test Clapeyron.lb_volume(model,T0,z3) ≈ 3*Clapeyron.lb_volume(model,T0,z)
    @test Clapeyron.lb_volume(model,T0,z7) ≈ 7*Clapeyron.lb_volume(model,T0,z)
    @test Clapeyron.T_scale(model,z3) ≈ Clapeyron.T_scale(model,z)
    @test Clapeyron.p_scale(model,z3) ≈ Clapeyron.p_scale(model,z)
end

function test_recombine(model,innermodels = nothing)
    model1 = deepcopy(model)
    model2 = deepcopy(model)
    Clapeyron.recombine!(model1)
    _test_recombine(model1,model2)
    if innermodels != nothing
        for field in innermodels
            innermodel = getfield(model1,field)
            innermodel2 = getfield(model2,field)
            _test_recombine(innermodel,innermodel2)
        end
    end
end

function _test_recombine(model1,model2)
    if hasfield(typeof(model1),:params)
        params1 = model1.params
        params2 = model2.params
        for i in 1:fieldcount(typeof(params1))
            p1 = getfield(params1,i)
            p2 = getfield(params2,i)
            if p1 isa SingleParam || p1 isa PairParam
                @test p1.values == p2.values
            elseif p1 isa AssocParam
                @test p1.values.values == p2.values.values
                @test p1.values.inner_indices == p2.values.inner_indices
                @test p1.values.outer_indices == p2.values.outer_indices
            elseif p1 isa Clapeyron.MixedGCSegmentParam
                @test p1.values.v == p2.values.v
            end
        end
    end
end

function test_kl(model; test_k = true, test_l = true)  
    model2 = deepcopy(model)
    
    if test_k
        k_orig = Clapeyron.get_k(model2)
        
        k0 = zeros(size(k_orig))
        Clapeyron.set_k!(model2,k0)
        @test k0 ≈ Clapeyron.get_k(model2)
        
        k1 = zeros(size(k_orig))
        s1,s2 = size(k1)
        for i in 1:s1
            for j in (i + 1):s2
                k1[j,i] = 0.01
            end
        end
        Clapeyron.set_k!(model2,k1)
        @test k1 ≈ Clapeyron.get_k(model2)

        Clapeyron.set_k!(model2,k_orig)
        @test k_orig ≈ Clapeyron.get_k(model2)
    end

    if test_l
        l_orig = Clapeyron.get_l(model2)
        
        l0 = zeros(size(l_orig))
        Clapeyron.set_l!(model2,l0)
        @test l0 ≈ Clapeyron.get_l(model2)
        
        l1 = zeros(size(l_orig))
        s1,s2 = size(l1)
        for i in 1:s1
            for j in (i + 1):s2
                l1[j,i] = 0.01
            end
        end
        Clapeyron.set_l!(model2,l1)
        @test l1 ≈ Clapeyron.get_l(model2)

        Clapeyron.set_l!(model2,l_orig)
        @test l_orig ≈ Clapeyron.get_l(model2)
    end
end
test_k(model) = test_kl(model,test_l = false)
test_l(model) = test_kl(model,test_k = false)

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
include_distributed("test_methods_eos.jl",5)
include_distributed("test_methods_api.jl",2)
include_distributed("test_methods_api_flash.jl",3)
include_distributed("test_methods_electrolytes.jl",1)
include_distributed("test_estimation.jl",1)
include_distributed("test_issues.jl",1)

