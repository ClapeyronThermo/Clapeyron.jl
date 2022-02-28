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

activity_coefficient(model::UNIFACModel,p,T,z) = activity_coefficient_ad(model,p,T,z)

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
            τji = exp(-evalpoly(T,(A[j,i],B[j,i],C[j,i]))*invT)
            ∑θpτ += θpj*τji
        end
        G_res += q_pi*X[i]*log(∑θpτ)
    end
    return -m̄*G_res
end

function excess_gibbs_free_energy(model::UNIFACModel,p,T,z)
    g_comp = excess_g_comb(model,p,T,z)
    g_res = excess_g_res(model,p,T,z)
    return (g_comp+g_res)*R̄*T 
end
