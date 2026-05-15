abstract type structSAFTgammaMieModel <: SAFTgammaMieModel end

struct structSAFTgammaMie{I,T} <: structSAFTgammaMieModel
    components::Vector{String}
    groups::GroupParam{T}
    sites::SiteParam
    params::SAFTgammaMieParam{T}
    idealmodel::I
    vrmodel::SAFTVRMie{I,T}
    epsilon_mixing::Symbol
    assoc_options::AssocOptions
    references::Array{String,1}
end

function structSAFTgammaMie(comps,groups,sites,params,idealmodel,pcsaftmodel,epsilon_mixing,assoc,refs)
    T = eltype(params)
    I = typeof(idealmodel)
    return structSAFTgammaMie{I,T}(comps,groups,sites,params,idealmodel,pcsaftmodel,epsilon_mixing,assoc,refs)
end

default_references(::Type{structSAFTgammaMie}) = ["10.1063/1.4851455", "10.1021/je500248h","10.1063/5.0048315", "doi.org/10.1021/acs.iecr.2c00198"]

"""
    structSAFTgammaMie <: SAFTgammaMieModel

    structSAFTgammaMie(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    epsilon_mixing = :default,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `vst`: Single Parameter (`Float64`) - Number of segments (no units)
- `S`: Single Parameter (`Float64`) - Shape factor for segment (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `shapefactor`: Single Parameter (`Float64`) - Shape factor for segment (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy `[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`
- `mixed_segment`: Mixed Group Contribution Parameter: ∑nᵢₖνₖmₖ

## Input models
- `idealmodel`: Ideal Model

## Description

s-SAFT-γ-Mie EoS

!!! note "Group Fragmentation"

    Molecule fragmentation into functional groups and connectivity information is available in GCIdentifier.jl, using `SAFTgammaMieGroups`.


## References
1. Shaahmadi,, F., Hurter, R.M., Burger, A.J., Cripwell, J.T. (2021). Improving the SAFT-γ Mie equation of state to account for functional group interactions in a structural (s-SAFT-γ Mie) framework: Linear and branched alkanes. The Journal of Chemical Physics, 154, 244102. [doi:10.1063/5.0048315 ](https://doi.org/10.1063/5.0048315 )
2. Schulze-Hulbe, A., Shaahmadi, F., Burger, A.J., Cripwell, J.T. (2022). Extending the Structural (s)-SAFT-γ Mie Equation of State to Primary Alcohols. Industrial & Engineering Chemistry Research, 61 (33), 12208-12228. [doi:10.1021/acs.iecr.2c00198](https://doi.org/10.1021/acs.iecr.2c00198)
"""
structSAFTgammaMie

function structSAFTgammaMie(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    group_userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    epsilon_mixing = :default,
    assoc_options = AssocOptions())

    epsilon_mixing = Symbol(epsilon_mixing)

    groups = GroupParam(components, ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv","SAFT/SAFTgammaMie/structSAFTgammaMie/structSAFTgammaMie_intragroups.csv"])
    params = getparams(groups, ["SAFT/SAFTgammaMie/structSAFTgammaMie","properties/molarmass_groups.csv"]; userlocations = userlocations, verbose = verbose)
    sites = params["sites"]
    components = groups.components
    
    segment = params["vst"]
    shapefactor = params["S"]
    mixed_segment = MixedGCSegmentParam(groups,shapefactor.values,segment.values)
    sigma = sigma_LorentzBerthelot(params["sigma"])
    sigma.values .*= 1E-10
    sigma3 = PairParam(sigma)
    sigma3 .= sigma3 .^ 3
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

    gcparams = SAFTgammaMieParam(segment,shapefactor,lambda_a,lambda_r,sigma,epsilon,epsilon_assoc,bondvol,mixed_segment)
    Mw_comps = group_sum(groups,params["Mw"])
    ideal_userlocations_updated = _update_idealuserlocations_for_GC(idealmodel,ideal_userlocations,Mw_comps)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations_updated,verbose)
    vrmodel = __SAFTVRMie(groups,gcparams,sites,idealmodel = init_idealmodel,assoc_options = assoc_options,epsilon_mixing = epsilon_mixing,verbose = verbose)
    vrmodel.params.Mw.values .= Mw_comps.values
    model = structSAFTgammaMie(components,groups,sites,gcparams,init_idealmodel,vrmodel,epsilon_mixing,assoc_options,default_references(structSAFTgammaMie))
    set_reference_state!(model,reference_state;verbose)
    return model
end

const sSAFTγMie = structSAFTgammaMie
@doc (@doc structSAFTgammaMie) sSAFTγMie
const sSAFTgammaMie = structSAFTgammaMie
@doc (@doc structSAFTgammaMie) sSAFTgammaMie

export structSAFTgammaMie,sSAFTgammaMie,sSAFTγMie

function a_chain(model::structSAFTgammaMieModel, V, T, z,_data = @f(data))
    _d,_,vrdata = _data
    _,ρS,_,_ζ_X,_ζst,_,_ = vrdata
    l = length(z)
    ∑z = ∑(z)

    ngroups = length(model.groups.flattenedgroups)

    _ϵ = model.params.epsilon
    _λr = model.params.lambda_r
    _λa = model.params.lambda_a
    _σ = model.params.sigma

    g_Mie = zero(V+T+first(z))*zeros(ngroups,ngroups)

    _KHS,ρS_∂KHS = @f(KHS_f_ρdf,_ζ_X)
    for k ∈ @groups
        ϵ = _ϵ[k,k]
        λa = _λa[k,k]
        λr = _λr[k,k] 
        σ = _σ[k,k]
        _C = @f(Cλ,λa,λr)
        dkk = _d[k]
        x_0ij = σ/dkk
        #calculations for a1 - diagonal
        aS_1_a,∂aS_1∂ρS_a = @f(aS_1_fdf,λa,_ζ_X,ρS)
        aS_1_r,∂aS_1∂ρS_r = @f(aS_1_fdf,λr,_ζ_X,ρS)
        B_a,∂B∂ρS_a = @f(B_fdf,λa,x_0ij,_ζ_X,ρS)
        B_r,∂B∂ρS_r = @f(B_fdf,λr,x_0ij,_ζ_X,ρS)

        #calculations for a2 - diagonal
        aS_1_2a,∂aS_1∂ρS_2a = @f(aS_1_fdf,2*λa,_ζ_X,ρS)
        aS_1_2r,∂aS_1∂ρS_2r = @f(aS_1_fdf,2*λr,_ζ_X,ρS)
        aS_1_ar,∂aS_1∂ρS_ar = @f(aS_1_fdf,λa+λr,_ζ_X,ρS)
        B_2a,∂B∂ρS_2a = @f(B_fdf,2*λa,x_0ij,_ζ_X,ρS)
        B_2r,∂B∂ρS_2r = @f(B_fdf,2*λr,x_0ij,_ζ_X,ρS)
        B_ar,∂B∂ρS_ar = @f(B_fdf,λr+λa,x_0ij,_ζ_X,ρS)
        α = _C*(1/(λa-3)-1/(λr-3))

        g_HSi = @f(g_HS,x_0ij,_ζ_X)
        #@show (g_HSi,i)
        ∂a_1∂ρ_S = _C*(x_0ij^λa*(∂aS_1∂ρS_a+∂B∂ρS_a)
                      - x_0ij^λr*(∂aS_1∂ρS_r+∂B∂ρS_r))
        #@show (∂a_1∂ρ_S,1)

        g_1_ = 3*∂a_1∂ρ_S-_C*(λa*x_0ij^λa*(aS_1_a+B_a)-λr*x_0ij^λr*(aS_1_r+B_r))
        #@show (g_1_,i)
        θ = exp(ϵ/T)-1
        γc = 10 * (-tanh(10*(0.57-α))+1) * _ζst*θ*exp(-6.7*_ζst-8*_ζst^2)
        ∂a_2∂ρ_S = 0.5*_C^2 *
            (ρS_∂KHS*(x_0ij^(2*λa)*(aS_1_2a+B_2a)
            - 2*x_0ij^(λa+λr)*(aS_1_ar+B_ar)
            + x_0ij^(2*λr)*(aS_1_2r+B_2r))
            + _KHS*(x_0ij^(2*λa)*(∂aS_1∂ρS_2a+∂B∂ρS_2a)
            - 2*x_0ij^(λa+λr)*(∂aS_1∂ρS_ar+∂B∂ρS_ar)
            + x_0ij^(2*λr)*(∂aS_1∂ρS_2r+∂B∂ρS_2r)))
    
        gMCA2 = 3*∂a_2∂ρ_S-_KHS*_C^2 *
        (λr*x_0ij^(2*λr)*(aS_1_2r+B_2r)-
            (λa+λr)*x_0ij^(λa+λr)*(aS_1_ar+B_ar)+
            λa*x_0ij^(2*λa)*(aS_1_2a+B_2a))
        g_2_ = (1+γc)*gMCA2
        #@show (g_2_,i)
        g_Mie[k,k] = g_HSi*exp(ϵ/T*g_1_/g_HSi+(ϵ/T)^2*g_2_/g_HSi)
        for l in 1:k-1
            ϵ = _ϵ[k,l]
            λa = _λa[k,l]
            λr = _λr[k,l] 
            σ = _σ[k,l]
            _C = @f(Cλ,λa,λr)
            dkl = (_d[k]+_d[l])/2
            x_0 = σ/dkl
            #calculations for a1 - diagonal
            aS_1_a,∂aS_1∂ρS_a = @f(aS_1_fdf,λa,_ζ_X,ρS)
            aS_1_r,∂aS_1∂ρS_r = @f(aS_1_fdf,λr,_ζ_X,ρS)
            B_a,∂B∂ρS_a = @f(B_fdf,λa,x_0,_ζ_X,ρS)
            B_r,∂B∂ρS_r = @f(B_fdf,λr,x_0,_ζ_X,ρS)

            #calculations for a2 - diagonal
            aS_1_2a,∂aS_1∂ρS_2a = @f(aS_1_fdf,2*λa,_ζ_X,ρS)
            aS_1_2r,∂aS_1∂ρS_2r = @f(aS_1_fdf,2*λr,_ζ_X,ρS)
            aS_1_ar,∂aS_1∂ρS_ar = @f(aS_1_fdf,λa+λr,_ζ_X,ρS)
            B_2a,∂B∂ρS_2a = @f(B_fdf,2*λa,x_0,_ζ_X,ρS)
            B_2r,∂B∂ρS_2r = @f(B_fdf,2*λr,x_0,_ζ_X,ρS)
            B_ar,∂B∂ρS_ar = @f(B_fdf,λr+λa,x_0,_ζ_X,ρS)
            α = _C*(1/(λa-3)-1/(λr-3))

            g_HSi = @f(g_HS,x_0,_ζ_X)
            #@show (g_HSi,i)
            ∂a_1∂ρ_S = _C*(x_0^λa*(∂aS_1∂ρS_a+∂B∂ρS_a)
                        - x_0^λr*(∂aS_1∂ρS_r+∂B∂ρS_r))
            #@show (∂a_1∂ρ_S,1)

            g_1_ = 3*∂a_1∂ρ_S-_C*(λa*x_0^λa*(aS_1_a+B_a)-λr*x_0^λr*(aS_1_r+B_r))
            #@show (g_1_,i)
            θ = exp(ϵ/T)-1
            γc = 10 * (-tanh(10*(0.57-α))+1) * _ζst*θ*exp(-6.7*_ζst-8*_ζst^2)
            ∂a_2∂ρ_S = 0.5*_C^2 *
                (ρS_∂KHS*(x_0^(2*λa)*(aS_1_2a+B_2a)
                - 2*x_0^(λa+λr)*(aS_1_ar+B_ar)
                + x_0^(2*λr)*(aS_1_2r+B_2r))
                + _KHS*(x_0^(2*λa)*(∂aS_1∂ρS_2a+∂B∂ρS_2a)
                - 2*x_0^(λa+λr)*(∂aS_1∂ρS_ar+∂B∂ρS_ar)
                + x_0^(2*λr)*(∂aS_1∂ρS_2r+∂B∂ρS_2r)))
        
            gMCA2 = 3*∂a_2∂ρ_S-_KHS*_C^2 *
            (λr*x_0^(2*λr)*(aS_1_2r+B_2r)-
                (λa+λr)*x_0^(λa+λr)*(aS_1_ar+B_ar)+
                λa*x_0^(2*λa)*(aS_1_2a+B_2a))
            g_2_ = (1+γc)*gMCA2
            #@show (g_2_,i)
            g_Mie[k,l] = g_HSi*exp(ϵ/T*g_1_/g_HSi+(ϵ/T)^2*g_2_/g_HSi)
            g_Mie[l,k] = g_Mie[k,l]
        end
    end

    S = model.params.shapefactor
    vˢ = model.params.segment
    v = model.groups.n_flattenedgroups
    v_intra = model.groups.n_intergroups

    achain = zero(V+T+first(z))

    for i ∈ @comps
        NBsi = sum(v[i].*vˢ)-1
        NBi = sum(v[i].*vˢ.*S)-1
        for k ∈ @groups
            bk = v[i][k]*(vˢ[k]-1)/NBsi
            achain+=z[i]*NBi*bk*log(g_Mie[k,k])
            for l ∈ 1:k
                bkl = v_intra[i][k,l]/NBsi
                achain+=z[i]*NBi*bkl*log(g_Mie[k,l])
            end
        end
    end

    return -achain/∑z
end