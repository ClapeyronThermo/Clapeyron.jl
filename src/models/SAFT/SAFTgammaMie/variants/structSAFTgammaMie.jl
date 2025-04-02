abstract type structSAFTgammaMieModel <: SAFTgammaMieModel end

struct structSAFTgammaMie{I,T} <: structSAFTgammaMieModel
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

s-SAFT-γ-Mie EoS

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

    groups = GroupParam(components, ["SAFT/SAFTgammaMie/SAFTgammaMie_groups.csv","SAFT/SAFTgammaMie/structSAFTgammaMie/structSAFTgammaMie_intragroups.csv"])
    params = getparams(groups, ["SAFT/SAFTgammaMie/structSAFTgammaMie","properties/molarmass_groups.csv"]; userlocations = userlocations, verbose = verbose)
    sites = params["sites"]
    components = groups.components
    
    gc_segment = params["vst"]
    shapefactor = params["S"]

    mw = group_sum(groups,params["Mw"])
    
    mixed_segment = MixedGCSegmentParam(groups,shapefactor.values,gc_segment.values)
    
    segment = SingleParam("segment",components,group_sum(groups,nothing))
    
    gc_sigma = sigma_LorentzBerthelot(params["sigma"])  
    gc_sigma.values .*= 1E-10
    gc_sigma3 = PairParam(gc_sigma)
    gc_sigma3.values .^= 3
    sigma3 = group_pairmean(mixed_segment,gc_sigma3)
    sigma3.values .= cbrt.(sigma3.values)
    sigma = sigma_LorentzBerthelot(sigma3)

    if epsilon_mixing == :default
        gc_epsilon = epsilon_HudsenMcCoubreysqrt(params["epsilon"], gc_sigma)
        epsilon = epsilon_HudsenMcCoubreysqrt(group_pairmean(mixed_segment,gc_epsilon),sigma)
    elseif epsilon_mixing == :hudsen_mccoubrey
        gc_epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], gc_sigma)
        epsilon = epsilon_HudsenMcCoubrey(group_pairmean(mixed_segment,gc_epsilon),sigma)
    else
        throw(error("invalid specification of ",error_color(epsilon_mixing),". available values are :default and :hudsen_mccoubrey"))
    end
    gc_lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    gc_lambda_r = lambda_LorentzBerthelot(params["lambda_r"])

    lambda_a = group_pairmean(mixed_segment,gc_lambda_a) |> lambda_LorentzBerthelot
    lambda_r = group_pairmean(mixed_segment,gc_lambda_r) |> lambda_LorentzBerthelot
 
    #GC to component model in association
    gc_epsilon_assoc = params["epsilon_assoc"]
    gc_bondvol = params["bondvol"]
    gc_bondvol,gc_epsilon_assoc = assoc_mix(gc_bondvol,gc_epsilon_assoc,gc_sigma,assoc_options,sites) #combining rules for association

    comp_sites = gc_to_comp_sites(sites,groups)
    comp_bondvol = gc_to_comp_sites(gc_bondvol,comp_sites)
    comp_epsilon_assoc = gc_to_comp_sites(gc_epsilon_assoc,comp_sites)

    gcparams = SAFTgammaMieParam(gc_segment, shapefactor,gc_lambda_a,gc_lambda_r,gc_sigma,gc_epsilon,gc_epsilon_assoc,gc_bondvol,mixed_segment)
    vrparams = SAFTVRMieParam(mw,segment,sigma,lambda_a,lambda_r,epsilon,comp_epsilon_assoc,comp_bondvol)
    
    idmodel = init_model(idealmodel,components,ideal_userlocations,verbose,reference_state)
    
    vr = SAFTVRMie(components,comp_sites,vrparams,idmodel,assoc_options,default_references(SAFTVRMie))
    γmierefs = ["10.1063/1.4851455", "10.1021/je500248h"]
    gmie = structSAFTgammaMie(components,groups,sites,gcparams,idmodel,vr,epsilon_mixing,assoc_options,γmierefs)
    set_reference_state!(gmie,reference_state;verbose)
    return gmie
end

const sSAFTγMie = structSAFTgammaMie
const sSAFTgammaMie = structSAFTgammaMie

export structSAFTgammaMie,sSAFTgammaMie,sSAFTγMie

function a_res(model::structSAFTgammaMieModel, V, T, z)
    _data = @f(data)
    dgc,X,vrdata = _data
    _,ρS,ζi,_ζ_X,_ζst,σ3x,m̄ = vrdata
    vrdata_disp = (dgc,ρS,ζi,_ζ_X,_ζst,σ3x,m̄)
    return @f(a_hs,_data) + a_disp(model,V,T,X,vrdata_disp)/sum(z) + @f(a_chain,_data) + @f(a_assoc,_data)
end

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