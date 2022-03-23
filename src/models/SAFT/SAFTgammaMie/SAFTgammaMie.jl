
include("utils.jl")
#just a holder for the z partitions.
#to allow split_model to work correctly

abstract type SAFTgammaMieModel <: SAFTVRMieModel end


struct SAFTgammaMieParam <: EoSParam
    segment::SingleParam{Int}
    shapefactor::SingleParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end


struct SAFTgammaMie{I,VR} <: SAFTgammaMieModel
    components::Vector{String}
    groups::GroupParam
    sites::SiteParam
    params::SAFTgammaMieParam
    idealmodel::I
    vrmodel::VR
    assoc_options::AssocOptions
    references::Array{String,1}
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

    _mw = group_sum(groups,params["Mw"])

    mix_segment!(groups,S,vst)
    segment = SingleParam("segment",components,group_sum(groups))

    gc_sigma = sigma_LorentzBerthelot(params["sigma"])
    gc_sigma.values .*= 1E-10
    gc_epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], gc_sigma)
    ϵ = gc_epsilon.values
    σ = gc_sigma.values
     #helper functions
     function ẑ(i, k)
        return v[i][k]*vst[k]*S[k] / ∑(v[i][l]*vst[l]*S[l] for l ∈ gc)
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
    gcparams = SAFTgammaMieParam(gc_segment, shapefactor,gc_lambda_a,gc_lambda_r,gc_sigma,gc_epsilon,gc_epsilon_assoc,gc_bondvol)
    vrparams = SAFTVRMieParam(segment,sigma,lambda_a,lambda_r,epsilon,comp_epsilon_assoc,comp_bondvol,_mw)
    idmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    vr = SAFTVRMie(vrparams, comp_sites, idmodel; ideal_userlocations, verbose, assoc_options)
    γmierefs = ["10.1063/1.4851455", "10.1021/je500248h"]
    gmie = SAFTgammaMie(components,groups,sites,gcparams,idmodel,vr,assoc_options,γmierefs)
    return gmie
end
@registermodel SAFTgammaMie

const SAFTγMie = SAFTgammaMie
export SAFTgammaMie,SAFTγMie

SAFTVRMie(model::SAFTgammaMie) = model.vrmodel

include("equations.jl")
