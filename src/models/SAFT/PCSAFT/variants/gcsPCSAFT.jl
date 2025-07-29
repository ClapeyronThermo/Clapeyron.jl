abstract type gcsPCSAFTModel <: sPCSAFTModel end

struct gcsPCSAFTParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    segment::SingleParam{T}
    msigma3::SingleParam{T}
    mepsilon::SingleParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T}
end

function gcsPCSAFTParam(Mw,m,mσ3,mϵ,ϵijab,β)
    return build_parametric_param(gcsPCSAFTParam,Mw,m,mσ3,mϵ,ϵijab,β)
end

struct gcsPCSAFT{I,T} <: gcsPCSAFTModel
    components::Vector{String}
    groups::GroupParam
    sites::SiteParam
    params::gcsPCSAFTParam{T}
    idealmodel::I
    pcsaftmodel::sPCSAFT{I,T}
    assoc_options::AssocOptions
    references::Array{String,1}
end

export gcsPCSAFT

"""
    gcsPCSAFT <: PCSAFTModel
    gcsPCSAFT(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `m`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) - Binary Interaction Paramater (no units). Interaction parameter is component-based.
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`
## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
## Input models
- `idealmodel`: Ideal Model
## Description
Group-contribution version of Simplified Perturbed-Chain SAFT (sPC-SAFT)
## References
1. Tihic, A., Kontogeorgis, G.M., von Solms, N., Michelsen, M.L. (2008). A predictive group-contribution simplified PC-SAFT equation of state: Application to polymer systems. Industrial & Engineering Chemistry Research, 47(15), 5092-5101. [doi:10.1021/ie0710768](https://doi.org/10.1021/ie0710768)
"""
gcsPCSAFT

function gcsPCSAFT(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

    groups = GroupParam(components,["SAFT/PCSAFT/gcsPCSAFT/gcsPCSAFT_groups.csv"]; group_userlocations = group_userlocations,verbose = verbose)
    params = getparams(groups, ["SAFT/PCSAFT/gcsPCSAFT/","properties/molarmass_groups.csv"]; userlocations = userlocations, verbose = verbose)
    components = groups.components
    k_params = getparams(components, ["SAFT/PCSAFT/sPCSAFT/sPCSAFT_unlike.csv"]; userlocations = userlocations, verbose = verbose)
    params["k"] = get(k_params,"k") do
        PairParam("k",components)
    end

    return gcsPCSAFT(groups,params;idealmodel,ideal_userlocations,reference_state,verbose,assoc_options)
end

default_references(::Type{gcsPCSAFT}) = ["10.1021/ie020753p"]

function gcsPCSAFT(groups::GroupParam, params::Dict{String,ClapeyronParam};
    idealmodel = BasicIdeal,
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

    components = groups.components

    sites = params["sites"]
    mw = params["Mw"]
    segment = params["segment"]
    msigma3 = params["msigma3"]
    msigma3.values .*= 1E-30
    mepsilon = params["mepsilon"]
    gc_epsilon_assoc = params["epsilon_assoc"]
    gc_bondvol = params["bondvol"]
    k = get(params,"k",nothing)
    
    bondvol,epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,nothing,assoc_options,sites) #combining rules for association
    gcparams = gcsPCSAFTParam(mw, segment, msigma3, mepsilon, epsilon_assoc, bondvol)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    pcmodel = sPCSAFT(groups,gcparams,sites;
                idealmodel = init_idealmodel,
                assoc_options = assoc_options,
                k = k,
                verbose = verbose)

    model = gcsPCSAFT(components,groups,sites,gcparams,pcmodel.idealmodel,pcmodel,assoc_options,default_references(gcsPCSAFT))
    set_reference_state!(model,reference_state;verbose)
    return model
end


function sPCSAFT(groups::GroupParam,
                param::gcsPCSAFTParam,
                sites::SiteParam = SiteParam(groups.flattenedgroups);
                idealmodel = BasicIdeal(),
                assoc_options = AssocOptions(),
                k = nothing,
                verbose = false)

    mw = group_sum(groups,param.Mw)
    segment = group_sum(groups,param.segment)

    sigma = group_sum(groups,param.msigma3)
    sigma.values ./= segment.values
    sigma.values .= cbrt.(sigma.values)
    sigma = sigma_LorentzBerthelot(sigma)

    epsilon = group_sum(groups,param.mepsilon)
    epsilon.values ./= segment.values
    epsilon = epsilon_LorentzBerthelot(epsilon,k)

    comp_sites = gc_to_comp_sites(sites,groups)
    bondvol = gc_to_comp_sites(param.bondvol,comp_sites)
    epsilon_assoc = gc_to_comp_sites(param.epsilon_assoc,comp_sites)

    pcparams = PCSAFTParam(mw, segment, sigma, epsilon, epsilon_assoc, bondvol)
    return sPCSAFT(groups.components, comp_sites, pcparams, idealmodel, assoc_options, default_references(sPCSAFT))
end

function recombine_impl!(model::gcsPCSAFTModel)
    groups = model.groups
    components = model.components
    sites = model.sites
    assoc_options = model.assoc_options
    pcparams = model.pcsaftmodel.params
    params = model.params

    segment = group_sum!(pcparams.segment,groups,params.segment)

    sigma = group_sum!(pcparams.sigma,groups,params.msigma3)
    diagvalues(sigma.values) ./= segment.values
    sigma.values .= cbrt.(sigma.values)
    sigma_LorentzBerthelot!(sigma)

    k = get_k(model.pcsaftmodel)
    epsilon = group_sum!(pcparams.epsilon,groups,params.mepsilon)
    diagvalues(epsilon.values) ./= segment.values
    epsilon_LorentzBerthelot!(epsilon,k)

    gc_bondvol,gc_epsilon_assoc = assoc_mix(params.bondvol,params.epsilon_assoc,nothing,assoc_options)
    params.bondvol.values.values[:] = gc_bondvol.values.values
    params.epsilon_assoc.values.values[:] = gc_epsilon_assoc.values.values

    comp_sites = gc_to_comp_sites(sites,groups)
    comp_bondvol = gc_to_comp_sites(gc_bondvol,comp_sites)
    comp_epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)

    pcparams.bondvol.values.values[:] = comp_bondvol.values.values
    pcparams.epsilon_assoc.values.values[:] = comp_epsilon_assoc.values.values
    return model
end

function a_res(model::gcsPCSAFTModel,V,T,z)
    return a_res(model.pcsaftmodel,V,T,z)
end

assoc_shape(model::gcsPCSAFTModel) = assoc_shape(model.pcsaftmodel)
getsites(model::gcsPCSAFTModel) = getsites(model.pcsaftmodel)

function lb_volume(model::gcsPCSAFTModel, z)
    return lb_volume(model.pcsaftmodel, z)
end

function T_scale(model::gcsPCSAFTModel, z)
    return T_scale(model.pcsaftmodel,z)
end

function T_scales(model::gcsPCSAFTModel)
    return T_scales(model.pcsaftmodel)
end

function p_scale(model::gcsPCSAFTModel,z)
    return p_scale(model.pcsaftmodel,z)
end