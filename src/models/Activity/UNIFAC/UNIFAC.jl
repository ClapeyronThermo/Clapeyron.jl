struct UNIFACParam <: EoSParam
    A::PairParam{Float64}
    B::PairParam{Float64}
    C::PairParam{Float64}
    R::SingleParam{Float64}
    Q::SingleParam{Float64}
    mixedsegment::SingleParam{Vector{Float64}}
end

abstract type UNIFACModel <: ActivityModel end

struct UNIFAC{c<:EoSModel} <: UNIFACModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    groups::GroupParam
    params::UNIFACParam
    puremodel::Vector{c}
    references::Array{String,1}
    unifac_cache::UNIFACCache
end

@registermodel UNIFAC
const modUNIFAC = UNIFAC
export UNIFAC

function UNIFAC(components; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    groups = GroupParam(components, ["Activity/UNIFAC/UNIFAC_groups.csv"]; verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/UNIFAC_like.csv", "Activity/UNIFAC/UNIFAC_unlike.csv"]; userlocations=userlocations, asymmetricparams=["A","B","C"], ignore_missing_singleparams=["A","B","C"], verbose=verbose)
    A  = params["A"]
    B  = params["B"]
    C  = params["C"]
    R  = params["R"]
    Q  = params["Q"]
    icomponents = 1:length(components)
    gc_mixedsegment = mix_segment(groups) #this function is used in SAFTγMie
    init_puremodel = [puremodel([groups.components[i]]) for i in icomponents]
    packagedparams = UNIFACParam(A,B,C,R,Q,gc_mixedsegment)
    references = String[]
    cache = UNIFACCache(groups,packagedparams)
    model = UNIFAC(components,icomponents,groups,packagedparams,init_puremodel,references,cache)
    return model
end

function activity_coefficient(model::UNIFACModel,V,T,z)
    return exp.(@f(lnγ_comb)+ @f(lnγ_res))
end

function lnγ_comb(model::UNIFACModel,V,T,z)
    x = z ./ sum(z)
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    q_p = model.unifac_cache.q_p
    Φ = r/dot(x,r)
    Φ_p = q_p/dot(x,q_p)
    θ = q/dot(x,q)
    lnγ_comb = @. log(Φ_p)+(1-Φ_p)-5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_comb
end

function lnγ_SG(model::UNIFACModel,V,T,z)

    x = z ./ sum(z)

    r =model.unifac_cache.r
    q =model.unifac_cache.q

    Φ = r/dot(x,r)
    θ = q/dot(x,q)
    lnγ_SG = @. -5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_SG
end

function lnγ_res(model::UNIFACModel,V,T,z)
    v  = model.groups.n_flattenedgroups
    _ψ = @f(Ψ)
    lnΓ_ = @f(lnΓ,_ψ)
    lnΓi_ = @f(lnΓi,_ψ)
    lnγ_res_ =  [sum(v[i][k].*(lnΓ_[k].-lnΓi_[i][k]) for k ∈ @groups) for i ∈ @comps]
    return lnγ_res_
end

function lnΓ(model::UNIFACModel,V,T,z,ψ = @f(ψ))
    Q = model.params.Q.values
    v  = model.groups.n_flattenedgroups
    x = z ./ sum(z)
    X = sum(v[i][:]*x[i] for i ∈ @comps) ./ sum(sum(v[i][k]*x[i] for k ∈ @groups) for i ∈ @comps)
    θ = X.*Q / dot(X,Q)
    lnΓ_ = Q.*(1 .-log.(sum(θ[m]*ψ[m,:] for m ∈ @groups)) .- sum(θ[m]*ψ[:,m]./sum(θ[n]*ψ[n,m] for n ∈ @groups) for m ∈ @groups))
    return lnΓ_
end

function lnΓi(model::UNIFACModel,V,T,z,ψ = @f(ψ))
    Q = model.params.Q.values
    v  = model.groups.n_flattenedgroups
    ψ = @f(Ψ)
    X = [v[i][:] ./ sum(v[i][k] for k ∈ @groups) for i ∈ @comps]
    θ = [X[i][:].*Q ./ sum(X[i][n]*Q[n] for n ∈ @groups) for i ∈ @comps]
    lnΓi_ = [Q.*(1 .-log.(sum(θ[i][m]*ψ[m,:] for m ∈ @groups)) .- sum(θ[i][m]*ψ[:,m]./sum(θ[i][n]*ψ[n,m] for n ∈ @groups) for m ∈ @groups)) for i ∈ @comps]
    return lnΓi_
end

function Ψ(model::UNIFACModel,V,T,z)
    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    return @. exp(-(A+B*T+C*T^2)/T)
end

function excess_g_SG(model::UNIFACModel,p,T,z=SA[1.0])
    _0 = zero(eltype(z))
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    n = sum(z)
    invn = 1/n
    Φm = dot(r,z)*invn
    θm =  dot(q,z)*invn
    G_comb = _0
    for i ∈ @comps
        xi = z[i]*invn
        Φi = r[i]/Φm
        θi = q[i]/θm
        G_comb += 5*q[i]*xi*log(θi/Φi)
    end
    return n*G_comb
end

function excess_g_comb(model::UNIFACModel,p,T,z=SA[1.0])
    _0 = zero(eltype(z))
    r =model.unifac_cache.r
    q =model.unifac_cache.q
    q_p = model.unifac_cache.q_p
    n = sum(z)
    invn = 1/n
    Φm = dot(r,z)*invn
    Φpm = dot(q_p,z)*invn
    θm =  dot(q,z)*invn
    G_comb = _0
    for i ∈ @comps
        xi = z[i]*invn
        Φi = r[i]/Φm
        Φpi = q_p[i]/Φpm
        θi = q[i]/θm
        G_comb += xi*log(Φpi) + 5*q[i]*xi*log(θi/Φi)
    end
    return n*G_comb
end

function excess_g_res(model::UNIFACModel,p,T,z=SA[1.0])
    _0 = zero(T+first(z))
    Q = model.params.Q.values
    A = model.params.A.values
    B = model.params.B.values
    C = model.params.C.values
    n = sum(z)
    invT = 1/T
    mi  = model.params.mixedsegment.values
    m = model.unifac_cache.m
    m̄ = dot(z,m)
    m̄inv = 1/m̄
    θpm = zero(eltype(z))
    for i in @groups
        θpm += dot(z,mi[i])*m̄inv*Q[i]
    end
    G_res = _0
    for i ∈ @groups
        q_pi = Q[i]
        xi = dot(z,mi[i])*m̄inv  
        ∑θpτ = _0
        for j ∈ @groups
            xj = dot(z,mi[j])*m̄inv    
            θpj = Q[j]*xj/θpm
            τji = exp(-evalpoly(T,(A[j,i],B[j,i],C[j,i]))*invT)
            ∑θpτ += θpj*τji
        end
        G_res += q_pi*xi*log(∑θpτ)*m̄
    end
    return -n*G_res
end

function excess_gibbs_free_energy(model::UNIFACModel,p,T,z)
    g_comp = excess_g_comb(model,p,T,z)
    g_res = excess_g_res(model,p,T,z)
    return (g_comp+g_res)*R̄*T 
end
