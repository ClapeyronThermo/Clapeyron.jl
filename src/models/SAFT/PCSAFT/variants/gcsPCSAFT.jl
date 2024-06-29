abstract type gcsPCSAFTModel <: sPCSAFTModel end

struct gcsPCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    msigma3::SingleParam{Float64}
    mepsilon::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

struct gcsPCSAFT{I,PC} <: gcsPCSAFTModel
    components::Vector{String}
    groups::GroupParam
    sites::SiteParam
    params::gcsPCSAFTParam
    idealmodel::I
    pcsaftmodel::PC
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
- `k`: Pair Parameter (`Float64`) - Binary Interaction Paramater (no units)
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
    gc_params = getparams(groups, ["SAFT/PCSAFT/gcsPCSAFT/","properties/molarmass_groups.csv"]; userlocations = userlocations, verbose = verbose)
    components = groups.components
    params = getparams(components, ["SAFT/PCSAFT/sPCSAFT/sPCSAFT_unlike.csv"]; userlocations = userlocations, verbose = verbose)
    
    sites = gc_params["sites"]

    gc_mw = gc_params["Mw"]
    mw = group_sum(groups,gc_mw)

    gc_segment = gc_params["segment"]
    segment = group_sum(groups,gc_segment)

    gc_msigma3 = gc_params["msigma3"]
    gc_msigma3.values .*= 1E-30
    sigma = group_sum(groups,gc_msigma3)
    sigma.values ./= segment.values
    sigma.values .= cbrt.(sigma.values)
    sigma = sigma_LorentzBerthelot(sigma)

    gc_mepsilon = gc_params["mepsilon"]
    epsilon = group_sum(groups,gc_mepsilon)
    epsilon.values ./= segment.values

    k = params["k"]
    epsilon = epsilon_LorentzBerthelot(epsilon,k)

    gc_epsilon_assoc = gc_params["epsilon_assoc"]
    gc_bondvol = gc_params["bondvol"]
    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,[],assoc_options,sites) #combining rules for association

    comp_sites = gc_to_comp_sites(sites,groups)
    bondvol = gc_to_comp_sites(gc_bondvol,comp_sites)
    epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)

    gcparams = gcsPCSAFTParam(gc_mw, gc_segment, gc_msigma3, gc_mepsilon, gc_epsilon_assoc,gc_bondvol)
    params = PCSAFTParam(mw, segment, sigma, epsilon, epsilon_assoc, bondvol)
    
    idmodel = init_model(idealmodel,components,ideal_userlocations,verbose,reference_state)

    references = ["10.1021/ie020753p"]
    pc = sPCSAFT(components,sites,params,idmodel, assoc_options, references)
    model = gcsPCSAFT(components, groups, sites, gcparams,idmodel,pc, assoc_options, references)
    return model
end

function recombine_impl!(model::gcsPCSAFTModel)
    groups = model.groups
    components = model.components
    sites = model.sites
    assoc_options = model.assoc_options

    gc_msigma3 = model.params.msigma3
    gc_mepsilon = model.params.mepsilon
    gc_segment = model.params.segment
    gc_epsilon_assoc = model.params.epsilon_assoc
    gc_bondvol = model.params.bondvol

    segment = group_sum(groups,gc_segment)
    model.pcsaftmodel.params.segment.values[:] = segment.values

    sigma = group_sum(groups,gc_msigma3)
    sigma.values ./= segment.values
    sigma.values .= cbrt.(sigma.values)
    sigma = sigma_LorentzBerthelot(sigma)
    model.pcsaftmodel.params.sigma.values[:] = sigma.values

    epsilon = group_sum(groups,gc_mepsilon)
    epsilon.values ./= segment.values
    epsilon = epsilon_LorentzBerthelot(epsilon)
    model.pcsaftmodel.params.epsilon.values[:] = epsilon.values

    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,gc_sigma,assoc_options)
    model.params.bondvol.values.values[:] = gc_bondvol.values.values
    model.params.epsilon_assoc.values.values[:] = gc_epsilon_assoc.values.values

    comp_sites = gc_to_comp_sites(sites,groups)
    comp_bondvol = gc_to_comp_sites(gc_bondvol,comp_sites)
    comp_epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)

    model.pcsaftmodel.params.bondvol.values.values[:] = comp_bondvol.values.values
    model.pcsaftmodel.params.epsilon_assoc.values.values[:] = comp_epsilon_assoc.values.values
    return model
end

function a_res(model::gcsPCSAFTModel,V,T,z)
    return a_res(model.pcsaftmodel,V,T,z)
end

assoc_shape(model::gcsPCSAFTModel) = assoc_shape(model.pcsaftmodel)
getsites(model::gcsPCSAFTModel) = getsites(model.pcsaftmodel)

function lb_volume(model::gcsPCSAFTModel, z = SA[1.0])
    return lb_volume(model.pcsaftmodel, z)
end

function T_scale(model::gcsPCSAFTModel,z=SA[1.0])
    return T_scale(model.pcsaftmodel,z)
end

function T_scales(model::gcsPCSAFTModel)
    return T_scales(model.pcsaftmodel)
end

function p_scale(model::gcsPCSAFTModel,z=SA[1.0])
    return p_scale(model.pcsaftmodel,z)
end