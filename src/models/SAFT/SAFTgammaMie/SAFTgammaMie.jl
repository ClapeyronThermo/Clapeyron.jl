
abstract type SAFTgammaMieModel <: SAFTVRMieModel end


struct SAFTgammaMieParam{T} <: ParametricEoSParam{T}
    segment::SingleParam{Int}
    shapefactor::SingleParam{T}
    lambda_a::PairParam{T}
    lambda_r::PairParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T}
    mixed_segment::MixedGCSegmentParam{T}
end

function SAFTgammaMieParam{T}(group::GroupParam,sites = nothing) where T <: Number
    gc = group.flattenedgroups
    segment = SingleParam{Int}("segment",gc)
    shapefactor = SingleParam{T}("shapefactor",gc)
    lambda_a = PairParam{T}("lambda_a",gc)
    lambda_r = PairParam{T}("lambda_r",gc)
    sigma = PairParam{T}("sigma",gc)
    epsilon = PairParam{T}("epsilon",gc)
    if sites isa SiteParam
        epsilon_assoc = AssocParam{T}("epsilon",sites)
        bondvol = AssocParam{T}("epsilon",sites)
    else
        epsilon_assoc = AssocParam{T}("epsilon",gc)
        bondvol = AssocParam{T}("epsilon",gc)
    end
    mixed_segment = MixedGCSegmentParam(group,shapefactor.values,segment.values)
    return SAFTgammaMieParam{T}(segment,shapefactor,lambda_a,lambda_r,sigma,epsilon,epsilon_assoc,bondvol,mixed_segment)
end

SAFTgammaMieParam(group::GroupParam,sites = nothing) = SAFTgammaMieParam{Float64}(group,sites)

function SAFTgammaMieParam(segment,shapefactor,lambda_a,lambda_r,sigma,epsilon,epsilon_assoc,bondvol,mixed_segment)
    t = (segment,shapefactor,lambda_a,lambda_r,sigma,epsilon,epsilon_assoc,bondvol,mixed_segment)
    return build_parametric_param(SAFTgammaMieParam,segment,shapefactor,lambda_a,lambda_r,sigma,epsilon,epsilon_assoc,bondvol,mixed_segment)
end

struct SAFTgammaMie{I,T} <: SAFTgammaMieModel
    components::Vector{String}
    groups::GroupParam
    sites::SiteParam
    params::SAFTgammaMieParam{T}
    idealmodel::I
    vrmodel::SAFTVRMie{I,T}
    epsilon_mixing::Symbol
    assoc_options::AssocOptions
    references::Array{String,1}
end



"""
    SAFTgammaMie <: SAFTModel

SAFTgammaMie(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    epsilon_mixing = :default,
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `vst`: Single Parameter (`Float64`) - Number of segments (no units)
- `S`: Single Parameter (`Float64`) - Shape factor for segment (no units)
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
- `mixed_segment`: Mixed Group Contribution Parameter: ∑nᵢₖνₖmₖ

## Input models
- `idealmodel`: Ideal Model

## Description

SAFT-γ-Mie EoS

!!! info
    You can choose between the Hudsen-McCoubrey combining rule (`√(ϵᵢ*ϵⱼ)*(σᵢ^3 * σⱼ^3)/σᵢⱼ^6`) or the default rule (`√(ϵᵢ*ϵⱼ*(σᵢ^3 * σⱼ^3))/σᵢⱼ^3`) by passing the `epsilon_mixing` argument.
    with arguments `:default` or `:hudsen_mccoubrey`

## References
1. Papaioannou, V., Lafitte, T., Avendaño, C., Adjiman, C. S., Jackson, G., Müller, E. A., & Galindo, A. (2014). Group contribution methodology based on the statistical associating fluid theory for heteronuclear molecules formed from Mie segments. The Journal of Chemical Physics, 140(5), 054107. [doi:10.1063/1.4851455](https://doi.org/10.1063/1.4851455)
2. Dufal, S., Papaioannou, V., Sadeqzadeh, M., Pogiatzis, T., Chremos, A., Adjiman, C. S., … Galindo, A. (2014). Prediction of thermodynamic properties and phase behavior of fluids and mixtures with the SAFT-γ Mie group-contribution equation of state. Journal of Chemical and Engineering Data, 59(10), 3272–3288. [doi:10.1021/je500248h](https://doi.org/10.1021/je500248h)

## List of available groups
|Name    |Description         |
|--------|--------------------|
|CH3     |Methyl              |
|CH2     |Methylene           |
|CH      |                    |
|C       |                    |
|aCH     |Aromatic CH         |
|aCCH2   |                    |
|aCCH    |                    |
|CH2=    |Alkene group        |
|CH=     |                    |
|cCH2    |Cyclic alkane group |
|COOH    |Carboxylic acid group|
|COO     |Ester group         |
|OH      |Hydroxyl            |
|CH2OH   |Methylene hydroxyl group|
|CHOH    |                    |
|NH2     |Amine group         |
|NH      |                    |
|N       |                    |
|cNH     |                    |
|cN      |                    |
|CH=     |                    |
|aCCH3   |                    |
|aCOH    |                    |
|cCH     |                    |
|cCHNH   |                    |
|cCHN    |                    |
|aCCOaC  |                    |
|aCCOOH  |                    |
|aCNHaC  |                    |
|CH3CO   |                    |
|eO      |End ether oxygen    |
|cO      |Center ether oxygen |
"""
SAFTgammaMie

function SAFTgammaMie(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    epsilon_mixing = :default,
    assoc_options = AssocOptions())

    groups = GroupParam(components, ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv"]; group_userlocations = group_userlocations,verbose = verbose)
    params = getparams(groups, ["SAFT/SAFTgammaMie","properties/molarmass_groups.csv"]; userlocations = userlocations, verbose = verbose)

    return SAFTgammaMie(groups, params;
                        idealmodel = idealmodel,
                        ideal_userlocations = ideal_userlocations,
                        reference_state = reference_state,
                        verbose = verbose,
                        epsilon_mixing = epsilon_mixing,
                        assoc_options = assoc_options)
end

function SAFTgammaMie(groups::GroupParam, params::Dict{String,ClapeyronParam};
    idealmodel = BasicIdeal,
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    epsilon_mixing = :default,
    assoc_options = AssocOptions())

    
    sites = params["sites"]
    components = groups.components
    
    segment = params["vst"]
    shapefactor = params["S"]
    mixed_segment = MixedGCSegmentParam{Base.eltype(shapefactor)}(groups,shapefactor.values,segment.values)
    sigma = sigma_LorentzBerthelot(params["sigma"])
    sigma.values .*= 1E-10
    sigma3 = PairParam(sigma)
    sigma3.values .^= 3
    lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    
    if epsilon_mixing == :default
        epsilon = epsilon_HudsenMcCoubreysqrt(params["epsilon"], sigma)
    elseif epsilon_mixing == :hudsen_mccoubrey
        epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], sigma)
    else
        throw(error("invalid specification of ",error_color(epsilon_mixing),". available values are :default and :hudsen_mccoubrey"))
    end

    #GC to component model in association
    bondvol0 = params["bondvol"]
    epsilon_assoc0 = params["epsilon_assoc"]
    bondvol,epsilon_assoc = assoc_mix(bondvol0,epsilon_assoc0,sigma,assoc_options,sites) #combining rules for association

    gcparams = SAFTgammaMieParam(segment, shapefactor,lambda_a,lambda_r,sigma,epsilon,epsilon_assoc,bondvol,mixed_segment)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    vrmodel = SAFTVRMie(groups,gcparams,sites,idealmodel = init_idealmodel,assoc_options = assoc_options,epsilon_mixing = epsilon_mixing,verbose = verbose)
    group_sum!(vrmodel.params.Mw,groups,params["Mw"])
    model = SAFTgammaMie(components,groups,sites,gcparams,init_idealmodel,vrmodel,epsilon_mixing,assoc_options,default_references(SAFTgammaMie))
    set_reference_state!(model,reference_state;verbose)
    return model
end

mw(model::SAFTgammaMieModel) = mw(model.vrmodel)
molecular_weight(model::SAFTgammaMieModel,z) = molecular_weight(model.vrmodel,z)

const SAFTγMie = SAFTgammaMie
export SAFTgammaMie,SAFTγMie

SAFTVRMie(model::SAFTgammaMieModel) = model.vrmodel

function SAFTVRMie(groups::GroupParam,param::SAFTgammaMieParam,sites::SiteParam = SiteParam(groups.flattenedgroups);
    idealmodel = BasicIdeal(),assoc_options = AssocOptions(),
    epsilon_mixing = :default,
    verbose = false)

    verbose && @info("SAFTγ-Mie: creating SAFTVRMie model from SAFTγ-Mie parameters.")
    components = groups.components
    mixed_segment = param.mixed_segment
    gc_segment = param.segment
    shapefactor = param.shapefactor
    gc_sigma = param.sigma
    gc_epsilon = param.epsilon
    gc_lambda_r = param.lambda_r
    gc_lambda_a = param.lambda_a
    gc_bondvol = param.bondvol
    gc_epsilon_assoc = param.epsilon_assoc

    #segment
    segment = SingleParam("segment",components,group_sum(mixed_segment,nothing))

    #sigma
    gc_sigma3 = PairParam(gc_sigma)
    gc_sigma3.values .^= 3
    sigma3 = group_pairmean(mixed_segment,gc_sigma3)
    sigma3.values .= cbrt.(sigma3.values)
    sigma = sigma_LorentzBerthelot(sigma3)

    #epsilon
    if epsilon_mixing == :default
        epsilon = epsilon_HudsenMcCoubreysqrt(group_pairmean(mixed_segment,gc_epsilon),sigma)
    elseif epsilon_mixing == :hudsen_mccoubrey
        epsilon = epsilon_HudsenMcCoubrey(group_pairmean(mixed_segment,gc_epsilon),sigma)
    else
        throw(error("invalid specification of ",error_color(epsilon_mixing),". available values are :default and :hudsen_mccoubrey"))
    end
    #lambda
    lambda_a = group_pairmean(mixed_segment,gc_lambda_a) |> lambda_LorentzBerthelot
    lambda_r = group_pairmean(mixed_segment,gc_lambda_r) |> lambda_LorentzBerthelot

    comp_sites = gc_to_comp_sites(sites,groups)
    comp_bondvol = gc_to_comp_sites(gc_bondvol,comp_sites)
    comp_epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)
    Mw = SingleParam("Mw",components, zeros(eltype(shapefactor), length(groups.components)))
    vrparams = SAFTVRMieParam(Mw,segment,sigma,lambda_a,lambda_r,epsilon,comp_epsilon_assoc,comp_bondvol)
    return SAFTVRMie(components,comp_sites,vrparams,idealmodel,assoc_options,default_references(SAFTVRMie))
end

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
    mixed_segment = model.params.mixed_segment

    mix_segment!(mixed_segment,groups,shapefactor.values,gc_segment.values)
    model.vrmodel.params.segment.values[:] = group_sum(mixed_segment,nothing)

    gc_sigma = sigma_LorentzBerthelot!(gc_sigma)
    gc_sigma3 = PairParam(gc_sigma)
    gc_sigma3.values .^= 3
    sigma3 = group_pairmean(mixed_segment,gc_sigma3)
    sigma3.values .= cbrt.(sigma3.values)
    comp_sigma = sigma_LorentzBerthelot(sigma3)
    model.vrmodel.params.sigma.values[:] = comp_sigma.values

    gc_epsilon = epsilon_HudsenMcCoubrey!(gc_epsilon, gc_sigma)
    if model.epsilon_mixing == :default
        gc_epsilon = epsilon_HudsenMcCoubreysqrt!(gc_epsilon, gc_sigma)
        comp_epsilon = epsilon_HudsenMcCoubreysqrt(group_pairmean(mixed_segment,gc_epsilon),model.vrmodel.params.sigma)
    else
        gc_epsilon = epsilon_HudsenMcCoubrey!(gc_epsilon, gc_sigma)
        comp_epsilon = epsilon_HudsenMcCoubrey(group_pairmean(mixed_segment,gc_epsilon),model.vrmodel.params.sigma)
    end
    model.vrmodel.params.epsilon.values[:] = comp_epsilon.values

    gc_lambda_a = lambda_LorentzBerthelot!(gc_lambda_a)
    gc_lambda_r = lambda_LorentzBerthelot!(gc_lambda_r)

    comp_lambda_a = group_pairmean(mixed_segment,gc_lambda_a) |> lambda_LorentzBerthelot
    model.vrmodel.params.lambda_a.values[:] = comp_lambda_a.values
    comp_lambda_r = group_pairmean(mixed_segment,gc_lambda_r) |> lambda_LorentzBerthelot
    model.vrmodel.params.lambda_r.values[:] = comp_lambda_r.values

    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,gc_sigma,assoc_options,sites)
    model.params.bondvol.values.values[:] = gc_bondvol.values.values
    model.params.epsilon_assoc.values.values[:] = gc_epsilon_assoc.values.values

    comp_sites = gc_to_comp_sites(sites,groups)
    comp_bondvol = gc_to_comp_sites(gc_bondvol,comp_sites)
    comp_epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)

    model.vrmodel.params.bondvol.values.values[:] = comp_bondvol.values.values
    model.vrmodel.params.epsilon_assoc.values.values[:] = comp_epsilon_assoc.values.values
    return model
end

default_references(::Type{SAFTgammaMie}) = ["10.1063/1.4851455", "10.1021/je500248h"]
