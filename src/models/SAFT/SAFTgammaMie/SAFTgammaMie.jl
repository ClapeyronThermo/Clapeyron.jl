
include("utils.jl")
#just a holder for the z partitions.
#to allow split_model to work correctly

abstract type SAFTgammaMieModel <: GCSAFTModel end

struct γMieZ <: EoSModel
    components::Vector{String}
    z::SingleParam{Vector{Float64}}
end

struct SAFTgammaMieParam <: EoSParam
    segment::SingleParam{Int}
    mixedsegment::SingleParam{Vector{Float64}}
    shapefactor::SingleParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

@registermodel γMieZ

struct SAFTgammaMie{I,VR} <: SAFTgammaMieModel
    components::Vector{String}
    groups::GroupParam
    sites::SiteParam
    params::SAFTgammaMieParam
    idealmodel::I
    vrmodel::VR
    mie_zfractions::γMieZ
    assoc_options::AssocOptions
    references::Array{String,1}
end

function split_model(model::SAFTgammaMie,subset = nothing)
    if !(subset === nothing)
        error("SAFTgammaMie does not support custom subsets")
    end
    pures = auto_split_model(model)
    for (i,pure) in pairs(pures)
        filter!(!iszero,pure.mie_zfractions.z.values[1])
        mixm = pure.params.mixedsegment.values
        for mi in mixm
            xi = mi[i]
            resize!(mi, 1)
            mi[1] = xi
        end
    end 
    return pures
end

function SAFTgammaMie(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())

    groups = GroupParam(components, ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; verbose=verbose)
    params,sites = getparams(groups, ["SAFT/SAFTgammaMie","properties/molarmass_groups.csv"]; userlocations=userlocations, verbose=verbose)
    _sites = sites.i_sites
    components = groups.components
    gc = groups.i_flattenedgroups
    comps = 1:length(components)
    gc_segment = params["vst"]
    shapefactor = params["S"]
    S = shapefactor.values
    vst = gc_segment.values
    v  = groups.n_flattenedgroups
    comp_segment = zeros(Float64,length(comps))
    for i ∈ comps
        res_i = 0.0
        vi = v[i]
        groups_i = groups.i_groups[i]
        for idx ∈ 1:length(groups_i)
            k = groups_i[idx]
            res_i += vi[k]*S[k]*vst[k]
        end
        comp_segment[i] = res_i
    end
    segment = SingleParam("segment",components,comp_segment)

    #used in x_S:
    #x_S(group i) = dot(z,mixsegment[i])/dot(z,m_vr)
    mixsegment =  [[v[i][k]*vst[k]*S[k] for i ∈ comps] for k ∈ gc]
    gc_mixedsegment = SingleParam("mixed segment",groups.flattenedgroups,mixsegment,[false for i ∈ gc],String[],String[])
    gc_sigma = sigma_LorentzBerthelot(params["sigma"])
    gc_sigma.values .*= 1E-10
    gc_epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], gc_sigma)
    ϵ = gc_epsilon.values
    σ = gc_sigma.values
     #helper functions
     function ẑ(i, k)
        return v[i][k]*vst[k]*S[k] / ∑(v[i][l]*vst[l]*S[l] for l ∈ gc)
    end

    zz = [[0.0 for k ∈ 1:length(v[i])] for i ∈ comps]
    zzparam = SingleParam("z fraction",components,zz,[false for i ∈ 1:length(components)],String[],String[])
    for i ∈ 1:length(zz)
        zzi = zz[i]
        for k ∈ 1:length(zzi)
            zzi[k] = ẑ(i, k)
        end
    end
    
    function σ̄(i)
        return cbrt(∑(∑(ẑ(i,k)*ẑ(i,l)*σ[k,l]^3 for l ∈ gc) for k ∈ gc))
    end
    
    function σ̄(i, j)
        return (σ̄(i) + σ̄(j))/2
    end

    function ϵ̄(i)
        return ∑(∑(ẑ(i,k)*ẑ(i,l)*ϵ[k,l] for l ∈ gc) for k ∈ gc)
    end
    
    function ϵ̄(i, j)
        if i == j
            return ϵ̄(i)
        else
            return sqrt(σ̄(i)*σ̄(j))/σ̄(i,j) * sqrt(ϵ̄(i)*ϵ̄(j))
        end
    end

    comp_ϵ = [ϵ̄(i, j) for (i,j) ∈ Iterators.product(comps,comps)]
    epsilon = PairParam("epsilon",components,comp_ϵ)
    
    comp_σ = [σ̄(i, j) for (i,j) ∈ Iterators.product(comps,comps)]
    sigma = PairParam("sigma",components,comp_σ)

    gc_lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    gc_lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    #@show gc_lambda_a
    gc_λa = gc_lambda_a.values
    gc_λr = gc_lambda_r.values
    
    function λi(ll,i)
        return ∑(∑(ẑ(i,k)*ẑ(i,l)*ll[k,l] for l ∈ gc) for k ∈ gc)
    end
    comp_lambda_a = [λi(gc_λa,i) for i ∈ comps]
    comp_lambda_r = [λi(gc_λr,i) for i ∈ comps]
    lambda_a = lambda_LorentzBerthelot(SingleParam("lambda_a",components,comp_lambda_a))
    lambda_r = lambda_LorentzBerthelot(SingleParam("lambda_r",components,comp_lambda_r))

    #assoc
    gc_epsilon_assoc = params["epsilon_assoc"]
    gc_bondvol = params["bondvol"]
    comp_sites,idx_dict = gc_to_comp_sites(sites,groups)
    assoc_idx = gc_to_comp_assoc_idx(gc_bondvol,comp_sites,idx_dict)
    assoc_idxs,outer,inner,outer_size,inner_size = assoc_idx.values,assoc_idx.outer_indices,assoc_idx.inner_indices,assoc_idx.outer_size,assoc_idx.inner_size
    _comp_bondvol = [gc_bondvol.values.values[i] for i ∈ assoc_idxs]
    _comp_epsilon_assoc = [gc_epsilon_assoc.values.values[i] for i ∈ assoc_idxs]
    compval_bondvol = Compressed4DMatrix(_comp_bondvol,outer,inner,outer_size,inner_size)
    compval_epsilon_assoc = Compressed4DMatrix(_comp_epsilon_assoc,outer,inner,outer_size,inner_size)
    comp_bondvol = AssocParam{Float64}("epsilon assoc",components,compval_bondvol,comp_sites.sites,String[],String[])
    comp_epsilon_assoc = AssocParam{Float64}("bondvol",components,compval_epsilon_assoc,comp_sites.sites,String[],String[])
    #mixing of Mw
    comp_mw = zeros(Float64,length(components))
    gc_mw = params["Mw"].values
    for i ∈ 1:length(components)
        vi = v[i]
        gi = groups.i_groups[i]
        mwi = 0.0
        for idx ∈ 1:length(gi)
            j = gi[idx]
            mwi += gc_mw[j]*vi[j]
        end
        comp_mw[i] =mwi
    end  
    _mw = SingleParam("molecular_weight",components,comp_mw)
    gcparams = SAFTgammaMieParam(gc_segment,gc_mixedsegment, shapefactor,gc_lambda_a,gc_lambda_r,gc_sigma,gc_epsilon,gc_epsilon_assoc,gc_bondvol)
    vrparams = SAFTVRMieParam(segment,sigma,lambda_a,lambda_r,epsilon,comp_epsilon_assoc,comp_bondvol,_mw)
    idmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    vr = SAFTVRMie(vrparams, comp_sites, idmodel; ideal_userlocations, verbose, assoc_options)
    mie_z = γMieZ(components,zzparam)
    γmierefs = ["10.1063/1.4851455", "10.1021/je500248h"]
    gmie = SAFTgammaMie(components,groups,sites,gcparams,idmodel,vr,mie_z,assoc_options,γmierefs)
    return gmie
end
@registermodel SAFTgammaMie

const SAFTγMie = SAFTgammaMie
export SAFTgammaMie,SAFTγMie

SAFTVRMie(model::SAFTgammaMie) = model.vrmodel

include("equations.jl")
