#=
utils.jl

each test file can be run by their own if you include this file first.
=#

using Test
using Clapeyron
using CoolProp
using Unitful
using Clapeyron.LinearAlgebra
using Clapeyron.StaticArrays
using Clapeyron.EstimationUtils
using Clapeyron: has_sites,has_groups
using Clapeyron: PH, QT, QP, TS, VT, PT
using Aqua

macro printline()  # useful in hunting for where tests get stuck
    file = split(string(__source__.file), Base.Filesystem.path_separator)[end]
    printstyled(">>", file, ":", __source__.line, "\n", color=:light_black)
end


function test_gibbs_duhem(model,V,T,z;rtol = 1e-14)
    for i in (2.0,3.0,5.0,7.0,11.0)
        a_res₀ = Clapeyron.a_res(model,V,T,z)
        @test a_res₀ ≈ Clapeyron.a_res(model,i*V,T,i*z) rtol = rtol
    end
    pures = Clapeyron.split_pure_model(model)
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
    eosmodel_is_approx(model1,model2)
    if innermodels != nothing
        for field in innermodels
            innermodel = getfield(model1,field)
            innermodel2 = getfield(model2,field)
            eosmodel_is_approx(innermodel,innermodel2)
        end
    end
end

function eosmodel_is_approx(model1,model2)
    if hasfield(typeof(model1),:params)
        params1 = model1.params
        params2 = model2.params
        for i in 1:fieldcount(typeof(params1))
            p1 = getfield(params1,i)
            p2 = getfield(params2,i)
            if p1 isa SingleParam || p1 isa PairParam
                if eltype(p1) <: Tuple 
                    @test stack(p1.values) ≈ stack(p2.values)
                else
                    @test p1.values ≈ p2.values
                end
            elseif p1 isa AssocParam
                @test p1.values.values ≈ p2.values.values
                @test p1.values.indices == p2.values.indices
                @test p1.values.site_offsets == p2.values.site_offsets
            elseif p1 isa Clapeyron.MixedGCSegmentParam
                @test p1.values.v ≈ p2.values.v
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

function test_repr(val;str = nothing,str_compact = nothing)
    io = IOBuffer()
    Base.show(io,MIME"text/plain"(),val)
    x = String(take!(io))
    @test !isempty(x)

    io_compact = IOBuffer()
    Base.show(io_compact,val)
    x_compact = String(take!(io_compact))
    @test !isempty(x_compact)
    if str != nothing
        for s in str
            @test occursin(s,x)
        end
    end

    if str_compact != nothing
        for s in str_compact
            @test occursin(s,x_compact)
        end
    end
end

function test_zero_alloc2(model)
    z = ones(length(model))/length(model)
    V = 0.03
    T = 300.0
    f(x) = Clapeyron.eos(x,V,T,z)
    f(model)
    @test @allocated(f(model)) == 0
end

function test_zero_alloc1(model)
    z = SA[1.0]
    V = 0.03
    T = 300.0
    f(x) = Clapeyron.eos(x,V,T,z)
    f(model)
    @test @allocated(f(model)) == 0
end

function test_assoc_matrix(K,tol = 1e-12)
    s1,s2 = size(K)
    @test s1 == s2
    x0 = ones(eltype(K),s1)
    
    x10 = Clapeyron.assoc_matrix_solve(K,x0)
    x1 = Clapeyron.assoc_matrix_solve(K)
    dx = x10 - x1
    @test sqrt(sum(abs2,dx)) < tol
end

function test_assoc_matrix(model,V,T,z,tol = 1e-12)
    K = Clapeyron.dense_assoc_site_matrix(model,V,T,z)
    test_assoc_matrix(K,tol)
end

#! format: off
function model500()
#created via eos_repr
return SAFTgammaMie{BasicIdeal, Float64}(
        ["DIMETHYLAMINE"], 
        GroupParam{Float64}(["DIMETHYLAMINE"], [["NH", "CH3"]], :unknown, [[1.0, 2.0]], [Matrix{Float64}(undef, 0, 0)], [[2, 1]], ["CH3", "NH"], [[2.0, 1.0]], String[]), 
        SiteParam(["CH3", "NH"], [String[], ["H", "e1"]], Clapeyron.PackedVofV([1, 1, 3], [1, 1]), [Int64[], [1, 2]], ["H", "e1", "e2"], [[0, 0, 0], [1, 1, 0]], [[0, 0, 0], [1, 2, 0]], String[], nothing), 
        Clapeyron.SAFTgammaMieParam{Float64}(
        SingleParam{Int64}("vst", ["CH3", "NH"], [1, 1], Bool[0, 0], ["/home/22796002/PhD/Smarts/SAFTgammaMie_like_Clapeyron.csv"], String[]), 
        SingleParam{Float64}("S", ["CH3", "NH"], [0.31696615176726495, 0.15785492572578758], Bool[0, 0], nothing, nothing), 
        PairParam{Float64}("lambda_a", ["CH3", "NH"], [6.0 6.0; 6.0 6.0], Bool[0 1; 1 0], nothing, nothing), 
        PairParam{Float64}("lambda_r", ["CH3", "NH"], [31.239630776055115 8.307492753603125; 8.307492753603125 18.449975142848835], Bool[0 0; 0 0], nothing, nothing), 
        PairParam{Float64}("sigma", ["CH3", "NH"], [5.245968421264455e-10 4.0342843067136425e-10; 4.0342843067136425e-10 2.822600192162831e-10], Bool[0 1; 1 0], nothing, nothing), 
        PairParam{Float64}("epsilon", ["CH3", "NH"], [467.0612223874854 268.98688435945206; 268.98688435945206 374.46738001552916], Bool[0 0; 0 0], nothing, nothing), 
        AssocParam{Float64}("epsilon_assoc", ["CH3", "NH"], 
        Clapeyron.Compressed4DMatrix([1693.6574275348516], [(2, 2)], [(1, 2)]), [String[], ["H", "e1"]], nothing, nothing), 
        AssocParam{Float64}("bondvol", ["CH3", "NH"], 
        Clapeyron.Compressed4DMatrix([4.532862604822759e-29], [(2, 2)], [(1, 2)]), [String[], ["H", "e1"]], nothing, nothing), 
        MixedGCSegmentParam{Float64}("mixed segment", ["DIMETHYLAMINE"], Clapeyron.PackedVofV([1, 3], [0.6339323035345299, 0.15785492572578758]))), 
        BasicIdeal(), 
        SAFTVRMie{BasicIdeal, Float64}(["DIMETHYLAMINE"], 
        SiteParam(["DIMETHYLAMINE"], [["NH/H", "NH/e1"]], Clapeyron.PackedVofV([1, 3], [1, 1]), [[21, 22]], ["CH3OH/H", "CH3OH/e1", "CO2/e1", "CO2/e2", "COOH/H", "COOH/e1", "COOH/e2", "CH2OH/H", "CH2OH/e1", "COO/e1", "OH/H", "OH/e1", "cNH/H", "cNH/e1", "NH2/H", "NH2/e1", "CH3CO/e1", "CHOH/H", "CHOH/e1", "N/e1", "NH/H", "NH/e1", "cN/e1", "aCCOOH/H", "aCCOOH/e1", "aCCOOH/e2", "aCH/e1", "aCOH/H", "aCOH/e1", "aCOH/e2", "aCCH3/e1", "aCCH2/e1", "H2O/H", "H2O/e1", "CH3COCH3/H", "CH3COCH3/e1", "CH3COCH3/e2"], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], String[], [[(23, 1), (23, 2)]]), 
        Clapeyron.SAFTVRMieParam{Float64}(
        SingleParam{Float64}("Mw", ["DIMETHYLAMINE"], [45.085], Bool[0], nothing, nothing), 
        SingleParam{Float64}("segment", ["DIMETHYLAMINE"], [0.7917872292603174], Bool[0], nothing, nothing), 
        PairParam{Float64}("sigma", ["DIMETHYLAMINE"], [4.854448918280263e-10;;], Bool[0;;], nothing, nothing), 
        PairParam{Float64}("lambda_a", ["DIMETHYLAMINE"], [6.000000000000001;;], Bool[0;;], nothing, nothing), 
        PairParam{Float64}("lambda_r", ["DIMETHYLAMINE"], [23.41048565515027;;], Bool[0;;], nothing, nothing), 
        PairParam{Float64}("epsilon", ["DIMETHYLAMINE"], [400.14816241181626;;], Bool[0;;], nothing, nothing), 
        AssocParam{Float64}("epsilon_assoc", ["DIMETHYLAMINE"], 
        Clapeyron.Compressed4DMatrix([1693.6574275348516], [(1, 1)], [(1, 2)]), [["NH/H", "NH/e1"]], nothing, nothing), 
        AssocParam{Float64}("bondvol", ["DIMETHYLAMINE"], 
        Clapeyron.Compressed4DMatrix([4.532862604822759e-29], [(1, 1)], [(1, 2)]), [["NH/H", "NH/e1"]], nothing, nothing)), 
        BasicIdeal(), 
        AssocOptions(1.0e-12, 1.0e-12, 1000, 0.5, :nocombining, false), ["10.1063/1.4819786", "10.1080/00268976.2015.1029027"]), :default, 
        AssocOptions(1.0e-12, 1.0e-12, 1000, 0.5, :nocombining, false), ["10.1063/1.4851455", "10.1021/je500248h"])
end
