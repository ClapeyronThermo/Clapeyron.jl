
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

"""
    SAFTVRSWModel <: SAFTModel

    SAFTVRSW(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    group_userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `shapefactor`: Single Parameter (`Float64`) - Shape factor for segment (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Model Parameters
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `shapefactor`: Single Parameter (`Float64`) - Shape factor for segment (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Input models
- `idealmodel`: Ideal Model

## Description

SAFT-γ-Mie EoS

## References
1. Papaioannou, V., Lafitte, T., Avendaño, C., Adjiman, C. S., Jackson, G., Müller, E. A., & Galindo, A. (2014). Group contribution methodology based on the statistical associating fluid theory for heteronuclear molecules formed from Mie segments. The Journal of Chemical Physics, 140(5), 054107. [doi:10.1063/1.4851455](https://doi.org/10.1063/1.4851455)
2. Dufal, S., Papaioannou, V., Sadeqzadeh, M., Pogiatzis, T., Chremos, A., Adjiman, C. S., … Galindo, A. (2014). Prediction of thermodynamic properties and phase behavior of fluids and mixtures with the SAFT-γ Mie group-contribution equation of state. Journal of Chemical and Engineering Data, 59(10), 3272–3288. [doi:10.1021/je500248h](https://doi.org/10.1021/je500248h)
"""
SAFTgammaMie

function SAFTgammaMie(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    group_userlocations = String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions(), kwargs...)

    groups = GroupParam(components, ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; group_userlocations = group_userlocations,verbose=verbose)
    params,sites = getparams(groups, ["SAFT/SAFTgammaMie","properties/molarmass_groups.csv"]; userlocations=userlocations, verbose=verbose)
    components = groups.components
    
    gc_segment = params["vst"]
    shapefactor = params["S"]

    mw = group_sum(groups,params["Mw"])
    
    mix_segment!(groups,shapefactor.values,gc_segment.values)
    
    segment = SingleParam("segment",components,group_sum(groups,nothing))
    
    gc_sigma = sigma_LorentzBerthelot(params["sigma"])  
    gc_sigma.values .*= 1E-10
    gc_sigma3 = PairParam(gc_sigma)
    gc_sigma3.values .^= 3
    sigma3 = group_pairmean(groups,gc_sigma3)
    sigma3.values .= cbrt.(sigma3.values)
    sigma = sigma_LorentzBerthelot(sigma3)

    gc_epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], gc_sigma)
    epsilon = epsilon_HudsenMcCoubrey(group_pairmean(groups,gc_epsilon),sigma)
    
    gc_lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    gc_lambda_r = lambda_LorentzBerthelot(params["lambda_r"])

    lambda_a = group_pairmean(groups,gc_lambda_a) |> lambda_LorentzBerthelot
    lambda_r = group_pairmean(groups,gc_lambda_r) |> lambda_LorentzBerthelot
 
    #GC to component model in association
    gc_epsilon_assoc = params["epsilon_assoc"]
    gc_bondvol = params["bondvol"]
    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,gc_sigma,assoc_options) #combining rules for association

    comp_sites,site_translator = gc_to_comp_sites(sites,groups)
    comp_bondvol = gc_to_comp_sites(gc_bondvol,comp_sites,site_translator)
    comp_epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites,site_translator)


    gcparams = SAFTgammaMieParam(gc_segment, shapefactor,gc_lambda_a,gc_lambda_r,gc_sigma,gc_epsilon,gc_epsilon_assoc,gc_bondvol)
    vrparams = SAFTVRMieParam(mw,segment,sigma,lambda_a,lambda_r,epsilon,comp_epsilon_assoc,comp_bondvol)
    
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

function recombine_impl!(model::SAFTgammaMieModel)
    groups = model.groups
    components = model.components
    sites = model.sites
    assoc_options = model.assoc_options

    gc_sigma = model.params.sigma
    gc_epsilon = model.params.epsilon
    gc_segment = model.params.segment
    shapefactor = model.params.shapefactor
    gc_lambda_r = model.params.lambda_r
    gc_lambda_a = model.params.lambda_a
    gc_epsilon_assoc = model.params.epsilon_assoc
    gc_bondvol = model.params.bondvol

    mix_segment!(groups,shapefactor.values,gc_segment.values)
    model.vrmodel.params.segment.values[:] = group_sum(groups,nothing)

    gc_sigma = sigma_LorentzBerthelot!(gc_sigma)
    gc_sigma3 = PairParam(gc_sigma)
    gc_sigma3.values .^= 3
    sigma3 = group_pairmean(groups,gc_sigma3)
    sigma3.values .= cbrt.(sigma3.values)
    comp_sigma = sigma_LorentzBerthelot(sigma3)
    model.vrmodel.params.sigma.values[:] = comp_sigma.values

    gc_epsilon = epsilon_HudsenMcCoubrey!(gc_epsilon, gc_sigma)
    comp_epsilon = epsilon_HudsenMcCoubrey(group_pairmean(groups,gc_epsilon),model.vrmodel.params.sigma)
    model.vrmodel.params.epsilon.values[:] = comp_epsilon.values

    gc_lambda_a = lambda_LorentzBerthelot!(gc_lambda_a)
    gc_lambda_r = lambda_LorentzBerthelot!(gc_lambda_r)

    comp_lambda_a = group_pairmean(groups,gc_lambda_a) |> lambda_LorentzBerthelot
    model.vrmodel.params.lambda_a.values[:] = comp_lambda_a.values
    comp_lambda_r = group_pairmean(groups,gc_lambda_r) |> lambda_LorentzBerthelot
    model.vrmodel.params.lambda_r.values[:] = comp_lambda_r.values

    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,gc_sigma,assoc_options)
    model.params.bondvol.values.values[:] = gc_bondvol.values.values
    model.params.epsilon_assoc.values.values[:] = gc_epsilon_assoc.values.values


    comp_sites,site_translator = gc_to_comp_sites(sites,groups)
    comp_bondvol = gc_to_comp_sites(gc_bondvol,comp_sites,site_translator)
    comp_epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites,site_translator)
    
    model.vrmodel.params.bondvol.values.values[:] = comp_bondvol.values.values
    model.vrmodel.params.epsilon_assoc.values.values[:] = comp_epsilon_assoc.values.values
    return model
end