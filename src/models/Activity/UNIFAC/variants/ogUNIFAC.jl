struct ogUNIFACParam <: EoSParam
    A::PairParam{Float64}
    R::SingleParam{Float64}
    Q::SingleParam{Float64}
    mixedsegment::SingleParam{Vector{Float64}}
end

abstract type ogUNIFACModel <: UNIFACModel end

struct ogUNIFAC{c<:EoSModel} <: ogUNIFACModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    groups::GroupParam
    params::ogUNIFACParam
    puremodel::Vector{c}
    references::Array{String,1}
    unifac_cache::UNIFACCache
end

@registermodel ogUNIFAC
export ogUNIFAC

function ogUNIFAC(components; puremodel=PR,
    userlocations=String[], 
     verbose=false)
    groups = GroupParam(components, ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_groups.csv"]; verbose=verbose)

    params = getparams(groups, ["Activity/UNIFAC/ogUNIFAC/ogUNIFAC_like.csv", "Activity/UNIFAC/ogUNIFAC/ogUNIFAC_unlike.csv"]; userlocations=userlocations, asymmetricparams=["A"], ignore_missing_singleparams=["A"], verbose=verbose)
    A  = params["A"]
    R  = params["R"]
    Q  = params["Q"]
    icomponents = 1:length(components)
    
    init_puremodel = [puremodel([groups.components[i]]) for i in icomponents]
    gc_mixedsegment = mix_segment(groups) #this function is used in SAFTγMie
    packagedparams = ogUNIFACParam(A,R,Q,gc_mixedsegment)
    references = String[]
    cache = UNIFACCache(groups,packagedparams)
    model = ogUNIFAC(components,icomponents,groups,packagedparams,init_puremodel,references,cache)
    return model
end

function lnγ_comb(model::ogUNIFACModel,V,T,z)
    Q = model.params.Q.values
    R = model.params.R.values

    v  = model.groups.n_flattenedgroups

    x = z ./ sum(z)

    r =[sum(v[i][k]*R[k] for k in @groups) for i in @comps]
    q =[sum(v[i][k]*Q[k] for k in @groups) for i in @comps]

    Φ = r/sum(x[i]*r[i] for i ∈ @comps)
    θ = q/sum(x[i]*q[i] for i ∈ @comps)
    lnγ_comb = @. log(Φ)+(1-Φ)-5*q*(log(Φ/θ)+(1-Φ/θ))
    return lnγ_comb
end

function Ψ(model::ogUNIFACModel,V,T,z)
    A = model.params.A.values
    return @. exp(-A/T)
end

function excess_g_comb(model::ogUNIFACModel,p,T,z=SA[1.0])
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
        G_comb += xi*log(Φi) + 5*q[i]*xi*log(θi/Φi)
    end
    return n*G_comb
end

function excess_g_res(model::ogUNIFACModel,p,T,z=SA[1.0])
    _0 = zero(T+first(z))
    Q = model.params.Q.values
    A = model.params.A.values
    invT = 1/T
    mi  = model.params.mixedsegment.values
    m̄ = dot(z,model.unifac_cache.m)
    m̄inv = 1/m̄
    X = [dot(z,mi_i)*m̄inv for mi_i in mi]
    θpm = dot(X,Q)
    G_res = _0
    for i ∈ @groups
        q_pi = Q[i]
        ∑θpτ = _0
        for j ∈ @groups
            θpj = Q[j]*X[j]/θpm
            τji = exp(-A[j,i]*invT)
            ∑θpτ += θpj*τji
        end
        G_res += q_pi*X[i]*log(∑θpτ)
    end
    return -m̄*G_res
end
