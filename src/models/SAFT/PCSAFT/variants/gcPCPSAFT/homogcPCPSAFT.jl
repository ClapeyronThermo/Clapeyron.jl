
abstract type HomogcPCPSAFTModel <: PCPSAFTModel end
struct HomogcPCPSAFT{I} <: HomogcPCPSAFTModel
    components::Vector{String}
    groups::GroupParam
    sites::SiteParam
    params::PCPSAFTParam
    idealmodel::I
    ppcmodel::PCPSAFT{I}
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
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `m`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `k`: Pair Parameter (`Float64`) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

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
    gc_params = getparams(groups,["SAFT/PCSAFT/gcPCPSAFT/homo/"]; userlocations = userlocations, verbose = verbose,ignore_missing_singleparams = ["k"])
    gc_components = components
    components = groups.components    
    sites = gc_params["sites"]
    
    gc_mw = gc_params["Mw"]
    mw = group_sum(groups,gc_mw)

    gc_segment = gc_params["segment"]
    segment = group_sum(groups,gc_segment)

    gc_sigma = gc_params["sigma"]
    gc_sigma.values .*= 1E-10
    gc_sigma.values .^= 3
    gc_sigma.values .*= gc_segment.values
    sigma = group_sum(groups,gc_sigma)
    sigma.values ./= segment.values
    sigma.values .= cbrt.(sigma.values)
    sigma = sigma_LorentzBerthelot(sigma)
    gc_k = get(gc_params, "k") do
        n_gc = length(gc_sigma.components)
        PairParam(gc_sigma.components, "k",zeros(n_gc,n_gc))
    end
    k = group_pairmean2(groups,gc_k)
    gc_epsilon = gc_params["epsilon"]
    gc_epsilon.values .*= gc_segment.values
    epsilon = group_sum(groups,gc_epsilon)
    epsilon.values ./= segment.values
    epsilon = epsilon_LorentzBerthelot(epsilon,k)

    gc_dipole = gc_params["dipole"]
    gc_dipole2 = SingleParam("Dipole squared", groups.flattenedgroups, gc_dipole.^2 ./ k_B*1e-36*(1e-10*1e-3))
    dipole2 = group_sum(groups,gc_dipole2)
    dipole2 = SingleParam("Dipole squared",components, dipole2 ./ segment)
    dipole = SingleParam("Dipole",components, sqrt.(dipole2 .* k_B ./ 1e-36 ./ (1e-10*1e-3)))

    gc_epsilon_assoc = gc_params["epsilon_assoc"]
    gc_bondvol = gc_params["bondvol"]

    comp_sites = gc_to_comp_sites(sites,groups)
    bondvol = gc_to_comp_sites(gc_bondvol,comp_sites)
    epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options,comp_sites) #combining rules for association

    gc_sigma = gc_params["sigma"] |> sigma_LorentzBerthelot
    gc_epsilon = gc_params["epsilon"] |> epsilon_LorentzBerthelot

    gcparams = PCPSAFTParam(gc_mw, gc_segment, gc_sigma, gc_epsilon, gc_dipole, gc_dipole2, gc_epsilon_assoc,gc_bondvol)
    params = PCPSAFTParam(mw, segment, sigma, epsilon, dipole, dipole2, epsilon_assoc, bondvol)
    
    idmodel = init_model(idealmodel,gc_components,ideal_userlocations,verbose,reference_state)

    references = ["10.1021/ie020753p"]
    pc = PCPSAFT(components,comp_sites,params,idmodel, assoc_options, references)
    model = HomogcPCPSAFT(components, groups, sites, gcparams,idmodel,pc, assoc_options, references)
    return model
end

function recombine_impl!(model::HomogcPCPSAFTModel)
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
    model.pcmodel.params.segment.values[:] = segment.values

    sigma = group_sum(groups,gc_msigma3)
    sigma.values ./= segment.values
    sigma.values .= cbrt.(sigma.values)
    sigma = sigma_LorentzBerthelot(sigma)
    model.pcmodel.params.sigma.values[:] = sigma.values

    epsilon = group_sum(groups,gc_mepsilon)
    epsilon.values ./= segment.values
    epsilon = epsilon_LorentzBerthelot(epsilon)
    model.pcmodel.params.epsilon.values[:] = epsilon.values

    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,gc_sigma,assoc_options)
    model.params.bondvol.values.values[:] = gc_bondvol.values.values
    model.params.epsilon_assoc.values.values[:] = gc_epsilon_assoc.values.values

    comp_sites = gc_to_comp_sites(sites,groups)
    comp_bondvol = gc_to_comp_sites(gc_bondvol,comp_sites)
    comp_epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)

    model.pcmodel.params.bondvol.values.values[:] = comp_bondvol.values.values
    model.pcmodel.params.epsilon_assoc.values.values[:] = comp_epsilon_assoc.values.values
    return model
end

function a_res(model::HomogcPCPSAFTModel,V,T,z)
    return a_res(model.ppcmodel,V,T,z)
end

assoc_shape(model::HomogcPCPSAFTModel) = assoc_shape(model.ppcmodel)
getsites(model::HomogcPCPSAFTModel) = getsites(model.ppcmodel)

function lb_volume(model::HomogcPCPSAFTModel, z = SA[1.0])
    return lb_volume(model.ppcmodel, z)
end

function T_scale(model::HomogcPCPSAFTModel,z=SA[1.0])
    return T_scale(model.ppcmodel,z)
end

function T_scales(model::HomogcPCPSAFTModel)
    return T_scales(model.ppcmodel)
end

function p_scale(model::HomogcPCPSAFTModel,z=SA[1.0])
    return p_scale(model.ppcmodel,z)
end

export HomogcPCPSAFT
