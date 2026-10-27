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

function fit_scaling(vals, is; exponents = (-3, -2, -1, 0, 1, 2, 3))
    v = vec(collect(vals))
    x = vec(collect(is))
    length(v) == length(x) ||
        throw(ArgumentError("vals and is must have the same length"))

    B = reduce(hcat, (x .^ m for m in exponents))   # n × k

    # Pivoted-QR least squares (this is what B \ v does for a tall B)
    c = B \ v                                        # k-vector

    # Rank by |c_i| · ||B_i|| (unchanged metric from your original)
    cn      = vec(sqrt.(sum(abs2, B, dims = 1)))     # k-vector
    contrib = abs.(c .* cn)                          # k-vector
    [(m = exponents[k], c = c[k]) for k in sortperm(contrib, rev = true)]
end

function test_gibbs_duhem(model, V, T, z; rtol = 1e-14, atol = 1e-15)
    is   = exp.(range(log(2), log(20); length=15))
    name = string(Clapeyron.parameterless_type(typeof(model)))

    a0   = Clapeyron.a_res(model, V, T, z)
    vals = [Clapeyron.a_res(model, i*V, T, i*z) for i in is]
    a_res_intensive   = all(abs.(vals .- a0) .≤ rtol * abs(a0) + atol) 
    a_res_intensive || @info "a_res scaling" fit_scaling(vals, is)
    @test a_res_intensive

    # SVector{1} (SingleComp) must give the same result as Vector{Float64}
    pures = Clapeyron.split_pure_model(model)
    singlecomp_dispatch_test = true
    bad = Tuple{Int64,Float64,Float64,Float64}[]
    for n in eachindex(pures), i in is
        a_sc = Clapeyron.a_res(pures[n], i*V, T, Clapeyron.SVector(i))
        a_v  = Clapeyron.a_res(pures[n], i*V, T, [i])
        abs(a_sc - a_v) ≤ rtol * max(abs(a_sc), abs(a_v)) + atol ||
            push!(bad, (n, i, a_sc, a_v))
    end
    singlecomp_dispatch_test = isempty(bad)
    singlecomp_dispatch_test || @info "SingleComp dispatch mismatch ($name)" bad
    @test singlecomp_dispatch_test
    !singlecomp_dispatch_test && (return nothing)

    pure_limit_test = true
    for n in eachindex(pures)
        x = zeros(length(model)); x[n] = 1.0
        for i in is
            a_mix  = Clapeyron.a_res(model,    i*V, T, i*x)
            a_pure = Clapeyron.a_res(pures[n], i*V, T, Clapeyron.SVector(i))
            abs(a_mix - a_pure) ≤ rtol * max(abs(a_mix), abs(a_pure)) + atol ||
                (pure_limit_test = false; break)
        end
    end
    pure_limit_test || @info "pure limit test for ($name) failed"
    @test pure_limit_test

end


function test_gibbs_duhem(model::Clapeyron.ActivityModel, V, T, z; rtol = 1e-14, atol = 1e-15)
    is   = exp.(range(log(2), log(20); length=15))
    name = string(Clapeyron.parameterless_type(typeof(model)))


    GE0   = Clapeyron.excess_gibbs_free_energy(model, V, T, z)
    Gvals = [Clapeyron.excess_gibbs_free_energy(model, i*V, T, i*z) for i in is]
    GEi   = is .* GE0
    GE_extensive    = all(abs.(Gvals .- GEi) .≤ rtol .* abs.(GEi) .+ atol)
    @test GE_extensive
    GE_extensive || @info "GE scaling" fit_scaling(Gvals, is)
    

    n   = length(model)
    bad = Tuple{Int64,Float64,Float64}[]
    for k in 1:n
        xk = zeros(n); xk[k] = 1.0
        for i in is
            g = Clapeyron.excess_gibbs_free_energy(model, i*V, T, i*xk)
            abs(g) ≤ atol || push!(bad, (k, i, g))
        end
    end
    pure_limit_GE_test = isempty(bad)
    @test pure_limit_GE_test
    pure_limit_GE_test || @info "nonzero pure-limit GE ($name)" bad
end

function test_volume(model,p,T,z = Clapeyron.SA[1.0];rtol = 1e-8,phase = :unknown)
    v = volume(model,p,T,z,phase = phase)
    @test p ≈ Clapeyron.pressure(model,v,T,z) rtol = rtol
    return v
end

function test_scales(model, T0 = 300.0; rtol = 1e-14, atol = 1e-15)
    is   = exp.(range(log(2), log(20); length=15))
    n    = length(model)
    z    = ones(n) ./ n
    name = Clapeyron.parameterless_type(typeof(model))

    v0   = Clapeyron.lb_volume(model, T0, z)
    vals = [Clapeyron.lb_volume(model, T0, i*z) for i in is]
    vi   = is .* v0
    lb_volume_consistency   = all(abs.(vals .- vi) .≤ rtol .* abs.(vi) .+ atol)
    @test lb_volume_consistency
    lb_volume_consistency || @info "$name lb_volume scaling" fit_scaling(vals, is)

    s0   = Clapeyron.T_scale(model, z)
    vals = [Clapeyron.T_scale(model, i*z) for i in is]
    T_scale_consistency   = all(abs.(vals .- s0) .≤ rtol * abs(s0) + atol)
    @test T_scale_consistency
    T_scale_consistency || @info "$name T_scale scaling" fit_scaling(vals, is)

    s0   = Clapeyron.p_scale(model, z)
    vals = [Clapeyron.p_scale(model, i*z) for i in is]
    p_scale_consistency   = all(abs.(vals .- s0) .≤ rtol * abs(s0) + atol)
    @test p_scale_consistency
    p_scale_consistency || @info "$name: p_scale scaling" fit_scaling(vals, is)
end

function test_recombine(model,innermodels = nothing)
    model1 = deepcopy(model)
    model2 = deepcopy(model)
    Clapeyron.recombine!(model1)
    eosmodel_is_approx(model1,model2)
    if innermodels !== nothing
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
    if str !== nothing
        for s in str
            @test occursin(s,x)
        end
    end

    if str_compact !== nothing
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

function test_excess_gibbs_free_energy(model::Clapeyron.ActivityModel,p,T,z)
    γ = Clapeyron.activity_coefficient(model,p,T,z)
    return Clapeyron.Rgas(model)*T*sum(z[i]*log(γ[i]) for i ∈ eachindex(z))
end

function test_activity_coefficient(model::Clapeyron.ActivityModel,p,T,z)
    X = Clapeyron.gradient_type(model,T+p,z)
    return exp.(Clapeyron.Solvers.gradient(x->Clapeyron.excess_gibbs_free_energy(model,p,T,x),z)/(Clapeyron.Rgas(model)*T))::X
end

struct VanLaar_GE <: Clapeyron.ActivityModel
    A12::Float64
    A21::Float64
end

function Clapeyron.excess_gibbs_free_energy(model::VanLaar_GE,p,T,z)
    A12 = model.A12 + 1e-5/T
    A21 = model.A21 + 2e-5/T
    ge = (A12*A21*z[1]*z[2]) / (A12*z[1] + A21*z[2])
    return Clapeyron.Rgas()*T*ge
end

struct VanLaar_lngamma <: Clapeyron.ActivityModel
    A12::Float64
    A21::Float64
end

function Clapeyron.lnγ_impl!(lnγ,model::VanLaar_lngamma,V,T,z)
    A12 = model.A12 + 1e-5/T
    A21 = model.A21 + 2e-5/T
    ax = A12*z[1] + A21*z[2]
    lnγ[1] = A12*(A21*z[2]/ax)^2
    lnγ[2] = A21*(A12*z[1]/ax)^2
    return lnγ
end

Base.length(::VanLaar_GE) = 2
Base.length(::VanLaar_lngamma) = 2

#used to reproduce old values of UNIFACFV and UNIFACFVPoly models
function activity_old_unifacfv(model::Clapeyron.UNIFACFVModel, p, T, z)
    RT = Rgas(model)*T
    lnγ_comb = Clapeyron.lnγ_comb_old(model, p, T, z)
    lnγ_res  = Clapeyron.ForwardDiff.gradient(n -> Clapeyron.excess_g_res(model, p, T, n), z) ./ RT
    lnγ_FV   = Clapeyron.ForwardDiff.gradient(n -> Clapeyron.excess_g_FV(model, p, T, n), z) ./ RT
    return exp.(lnγ_comb .+ lnγ_res .+ lnγ_FV)
end

function testutils_fd_derivative(f,i)
    h = sqrt(eps(i))
    return (1*f(i-2h)-8*f(i-1h)+0*f(i)+8*f(i+1h)-1*f(i+2h))/(12*h)
end

function testutils_fd_gradient(f, x::AbstractVector; method =:central, eps  =1e-6)
    n = length(x)
    g = similar(x, float(eltype(x)))
    if method == :central
        for i in 1:n
            xp = copy(x); xp[i] += eps
            xm = copy(x); xm[i] -= eps
            g[i] = (f(xp) - f(xm)) / (2eps)
        end
    elseif method == :forward
        f0 = f(x)
        for i in 1:n
            xp = copy(x); xp[i] += eps
            g[i] = (f(xp) - f0) / eps
        end
    elseif method == :backward
        f0 = f(x)
        for i in 1:n
            xm = copy(x); xm[i] -= eps
            g[i] = (f0 - f(xm)) / eps
        end
    else
        error("method must be :central, :forward, or :backward")
    end
    return g
end

function pure_generic_test(model)
    name = Clapeyron.model_showname(model)
    Tc,Pc,Vc = crit_pure(model)
    @test Tc > 0
    @test Pc > 0
    @test Vc > 0
    z1 = Clapeyron.SA[1.0]
    T1 = 0.9*Tc
    ps,vlsat,vvsat = saturation_pressure(model,T1)
    
    @test !isnan(ps) && !isnan(vlsat) && !isnan(vvsat)
    if isnan(ps) || isnan(vlsat) || isnan(vvsat) || isnan(T1)
        sat = (ps,vlsat,vvsat)
        @info "$name failed calculation of saturation pressure at T = $T1, got $sat"
    end
    vl1 = test_volume(model,ps,T1,z1,phase = :l)
    vv1 = test_volume(model,ps,T1,z1,phase = :v)
    @test vlsat ≈ vl1 rtol = 1e-6
    @test vvsat ≈ vv1 rtol = 1e-6
    Ts,_,_ = saturation_temperature(model,ps)
    @test !isnan(Ts)
    if isnan(Ts)
        @info "$name failed calculation of saturation temperature at p = $ps"
    end
    @test Ts ≈ T1 rtol = 1e-6
    B2 = second_virial_coefficient(model,Ts,z1)
    @test B2 < 0
    if isnan(B2)
        @info "$name got NaN second virial coefficient at T = $Ts"
    end
    
    is   = exp.(range(log(2), log(20); length=15))
    name = Clapeyron.parameterless_type(typeof(model))

    vals = [Clapeyron.second_virial_coefficient(model, Ts, i*z1) for i in is]
    vi   = is .* B2
    B2_consistency   = all(abs.(vals .- vi) .≤ 1e-6 .* abs.(vi) .+ 1e-10)
    @test B2_consistency

    if !B2_consistency
        if any(isnan,vals)
            @info "$name got NaN second virial coefficient at T = $Ts"
        else
            @info "$name second virial coefficient scaling" fit_scaling(vals, is)
        end
    end
end