
abstract type HomogcPCPSAFTModel <: PCPSAFTModel end

struct HomogcPCPSAFTParam{T} <: ParametricEoSParam{T}
    Mw::SingleParam{T}
    segment::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
    k::PairParam{T}
    dipole::SingleParam{T}
    dipole2::SingleParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T}
end

function HomogcPCPSAFTParam(Mw,segment,sigma,epsilon,k,dipole,dipole2,epsilon_assoc,bondvol)
    return build_parametric_param(HomogcPCPSAFTParam,Mw,segment,sigma,epsilon,k,dipole,dipole2,epsilon_assoc,bondvol)
end

struct HomogcPCPSAFT{I,T} <: HomogcPCPSAFTModel
    components::Vector{String}
    groups::GroupParam{T}
    sites::SiteParam
    params::HomogcPCPSAFTParam{T}
    idealmodel::I
    pcpmodel::PCPSAFT{I,T}
    assoc_options::AssocOptions
    references::Array{String,1}
end

export HomogcPCPSAFT

"""
    HomogcPCPSAFTModel <: PCPSAFTModel

    HomogcPCPSAFT(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `m`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `k`: Pair Parameter (`Float64`) - Binary Interaction Parameter (no units)
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `k`: Pair Parameter (`Float64`) - Binary Interaction Parameter (no units)
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Input models
- `idealmodel`: Ideal Model

## Description
Homosegmented Group-contribution Polar Perturbed-Chain SAFT (Homo-gc-PCP-SAFT)

## References
1. Sauer, E., Stavrou, M., Gross, J. (2014). Comparison between a Homo- and a Heterosegmented Group Contribution Approach Based on the Perturbed-Chain Polar Statistical Associating Fluid Theory Equation of State. Industrial & Engineering Chemistry Research, 53(38), 14854-14864. [doi:10.1021/ie502203w](https://doi.org/10.1021/ie502203w)
"""
HomogcPCPSAFT

function HomogcPCPSAFT(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    verbose = false,
    reference_state = nothing,
    assoc_options = AssocOptions(combining = :cr1))
    groups = GroupParam(components,["SAFT/PCSAFT/gcPCPSAFT/homo/HomogcPCPSAFT_groups.csv"]; group_userlocations = group_userlocations,verbose = verbose)
    params = getparams(groups,["SAFT/PCSAFT/gcPCPSAFT/homo/"]; userlocations = userlocations, verbose = verbose,ignore_missing_singleparams = ["k"])
    return HomogcPCPSAFT(groups,params;idealmodel,ideal_userlocations,reference_state,verbose,assoc_options)
end

default_references(::Type{HomogcPCPSAFT}) = references = ["10.1021/ie020753p"]

function HomogcPCPSAFT(groups::GroupParam,params::Dict{String,ClapeyronParam};
    idealmodel = BasicIdeal,
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

    components = groups.components
    sites = params["sites"]
    
    mw = params["Mw"]
    
    segment = params["segment"]
    
    _sigma = params["sigma"]
    _sigma.values .*= 1e-10
    
    sigma = sigma_LorentzBerthelot(_sigma)
    epsilon = params["epsilon"] |> epsilon_LorentzBerthelot
    
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    
    dipole = params["dipole"]
    dipole2 = SingleParam("Dipole squared", groups.flattenedgroups, dipole.^2 ./ k_B*1e-36*(1e-10*1e-3))

    k = get(params, "k") do
        n_gc = length(groups.n_flattenedgroups)
        PairParam(components,"k",zeros(n_gc,n_gc),fill(true,(n_gc,n_gc)))
    end

    gcparams = HomogcPCPSAFTParam(mw, segment, sigma, epsilon, k, dipole, dipole2, epsilon_assoc, bondvol)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    pcpmodel = PCPSAFT(groups,gcparams,sites;
                idealmodel = init_idealmodel,
                assoc_options = assoc_options,
                verbose = verbose)

    model = HomogcPCPSAFT(components,groups,sites,gcparams,pcpmodel.idealmodel,pcpmodel,assoc_options,default_references(HomogcPCPSAFT))
    set_reference_state!(model,reference_state;verbose)
    return model
end

function PCPSAFT(groups::GroupParam,
                param::HomogcPCPSAFTParam,
                sites::SiteParam = SiteParam(groups.flattenedgroups);
                idealmodel = BasicIdeal(),
                assoc_options = AssocOptions(),
                verbose = false)


    components = groups.components
    mw = group_sum(groups,param.Mw)

    #m(i) = ∑n(ik)*m(ik)
    segment = group_sum(groups,param.segment)
    
    #m(i)*(σ(i))^3 = ∑n(ik)*m(ik)*σ(ik)^3
    gc_msigma3 = param.sigma .^3 .* param.segment
    sigma = SingleParam("sigma",components,group_sum(groups,gc_msigma3))
    sigma.values ./= segment.values
    sigma.values .= cbrt.(sigma.values)
    sigma = sigma_LorentzBerthelot(sigma)

    #m(i)*ϵ(i) = ∑n(ik)*m(ik)*ϵ(ik)
    k = group_pairmean2(groups,param.k)
    gc_mepsilon = diagvalues(param.epsilon) .* param.segment.values
    epsilon = SingleParam("epsilon",components,group_sum(groups,gc_mepsilon))
    epsilon.values ./= segment.values
    epsilon = epsilon_LorentzBerthelot(epsilon,k)

    _dipole2 = group_sum(groups,param.dipole2.values)
    dipole2 = SingleParam("Dipole squared",components, _dipole2 ./ segment.values)
    dipole = SingleParam("Dipole",components, sqrt.(dipole2 .* k_B ./ 1e-36 ./ (1e-10*1e-3)))

    comp_sites = gc_to_comp_sites(sites,groups)
    bondvol = gc_to_comp_sites(param.bondvol,comp_sites)
    epsilon_assoc = gc_to_comp_sites(param.epsilon_assoc,comp_sites)

    pcpparams = PCPSAFTParam(mw, segment, sigma, epsilon, dipole, dipole2, epsilon_assoc, bondvol)

    return PCPSAFT(groups.components, comp_sites, pcpparams, idealmodel, assoc_options, default_references(PCPSAFT))
end


function recombine_impl!(model::HomogcPCPSAFTModel)
    groups = model.groups
    components = model.components
    sites = model.sites
    assoc_options = model.assoc_options
    gcparams = model.params
    pcpmodel = model.pcpmodel
    params = pcpmodel.params

    #recombine outer params
    sigma_LorentzBerthelot!(gcparams.sigma)
    epsilon_LorentzBerthelot!(gcparams.epsilon)

    gcparams.dipole2 .= gcparams.dipole.^2 ./ k_B*1e-36*(1e-10*1e-3)
    group_sum!(params.dipole2,groups,gcparams.dipole2)
    params.dipole2.values ./= params.segment 
    params.dipole .= sqrt.(params.dipole2 .* k_B ./ 1e-36 ./ (1e-10*1e-3))
 
    #recombine inner PCP model
    mw = group_sum!(params.Mw,groups,gcparams.Mw)
    segment = group_sum!(params.segment,groups,gcparams.segment)

    gc_msigma3 = gcparams.sigma .^3 .* gcparams.segment
    sigma_diag = diagvalues(params.sigma)
    group_sum!(sigma_diag,groups,gc_msigma3)
    sigma_diag ./= segment.values
    sigma_diag .= cbrt.(sigma_diag)
    sigma_LorentzBerthelot!(params.sigma)

    k = group_pairmean2(groups,gcparams.k)
    gc_mepsilon = diagvalues(gcparams.epsilon.values) .* gcparams.segment.values
    epsilon = group_sum!(params.epsilon,groups,gc_mepsilon)
    diagvalues(epsilon.values) ./= segment.values
    epsilon_LorentzBerthelot!(epsilon,k)

    group_sum!(params.dipole2,groups,gcparams.dipole2 ./ gcparams.segment)
    params.dipole .= sqrt.(params.dipole2 .* k_B ./ 1e-36 ./ (1e-10*1e-3))

    comp_sites = gc_to_comp_sites(sites,groups)
    comp_bondvol = gc_to_comp_sites(gcparams.bondvol,comp_sites)
    comp_epsilon_assoc = gc_to_comp_sites(gcparams.epsilon_assoc,comp_sites)

    bondvol,epsilon_assoc = assoc_mix(comp_bondvol,comp_epsilon_assoc,params.sigma,assoc_options,comp_sites)
    params.bondvol.values.values[:] = bondvol.values.values
    params.epsilon_assoc.values.values[:] = epsilon_assoc.values.values

    return model
end

function a_res(model::HomogcPCPSAFTModel,V,T,z)
    return a_res(model.pcpmodel,V,T,z)
end

assoc_shape(model::HomogcPCPSAFTModel) = assoc_shape(model.pcpmodel)
getsites(model::HomogcPCPSAFTModel) = getsites(model.pcpmodel)

function lb_volume(model::HomogcPCPSAFTModel, T, z)
    return lb_volume(model.pcpmodel, T, z)
end

function T_scale(model::HomogcPCPSAFTModel,z)
    return T_scale(model.pcpmodel,z)
end

function T_scales(model::HomogcPCPSAFTModel)
    return T_scales(model.pcpmodel)
end

function p_scale(model::HomogcPCPSAFTModel,z)
    return p_scale(model.pcpmodel,z)
end

export HomogcPCPSAFT
