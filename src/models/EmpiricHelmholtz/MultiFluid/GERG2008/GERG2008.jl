struct GERG2008 <: MultiFluidModel
    components::Vector{String}
    properties::MultiFluidPropertyParam
    ideal::MultiFluidIdealParam
    single::MultiFluidSingleParam
    pair::MultiFluidPairParam
    references::Vector{String}
end

@registermodel GERG2008

function GERG2008(components::Vector{String})
    I = getnames_gerg2008(components)
    params = getparams_gerg2008()
    Mw = each_split_model(params[:Mw],I)
    rhoc = each_split_model(params[:rhoc],I)
    vc = SingleParam("critical volume (dm3/mol)",rhoc.components,1 ./ rhoc.values)
    Tc = each_split_model(params[:Tc],I)
    pc = each_split_model(params[:pc],I)
    lb_v = each_split_model(params[:lb_v],I)
    properties = MultiFluidPropertyParam(Mw,rhoc,vc,Tc,pc,lb_v)
    n0 = each_split_model(params[:n0],I)
    theta = each_split_model(params[:theta],I)
    gamma_T = each_split_model(params[:gamma_T],I)
    gamma_v = each_split_model(params[:gamma_v],I)
    beta_T = each_split_model(params[:beta_T],I)
    beta_v = each_split_model(params[:beta_v],I)
    ideal = MultiFluidIdealParam(n0,theta,gamma_T,gamma_v,beta_T,beta_v)
    n = each_split_model(params[:n],I) |> pack_vectors
    t = each_split_model(params[:t],I) |> pack_vectors
    d = each_split_model(params[:d],I) |> pack_vectors
    c = each_split_model(params[:c],I) |> pack_vectors
    single = MultiFluidSingleParam(n,t,d,c)
    nij = each_split_model(params[:nij],I) |> pack_vectors
    tij = each_split_model(params[:tij],I) |> pack_vectors
    dij = each_split_model(params[:dij],I) |> pack_vectors
    Fij = each_split_model(params[:Fij],I)
    eta_ij = each_split_model(params[:eta_ij],I) |> pack_vectors
    epsilon_ij = each_split_model(params[:epsilon_ij],I) |> pack_vectors
    beta_ij = each_split_model(params[:beta_ij],I) |> pack_vectors
    gamma_ij = each_split_model(params[:gamma_ij],I) |> pack_vectors
    pair = MultiFluidPairParam(nij,tij,dij,Fij,beta_ij,gamma_ij,eta_ij,epsilon_ij)
    references = ["10.1021/je300655b"]
    return GERG2008(components,properties,ideal,single,pair,references)
end

export GERG2008

"""
    GERG2008 <: MultiFluidModel
    GERG2008(components::Vector{String})

## Imput Parameters

None

## Description

The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures. valid for 21 compounds (`Clapeyron.GERG2008_names`). 
```

a = a⁰ + aʳ

a⁰ = ∑xᵢ(a⁰ᵢ(τᵢ,δᵢ) + ln(xᵢ))
δᵢ = ρ/ρcᵢ
τᵢ = Tcᵢ/T
a⁰ᵢ = ln(δᵢ) + R∗/R[n⁰ᵢ₋₁ + n⁰ᵢ₋₂τᵢ + n⁰ᵢ₋₃ln(τᵢ) + ∑n⁰ᵢ₋ₖln(abs(sinh(ϑ₀ᵢ₋ₖτᵢ))) + ∑n⁰ᵢ₋ₖln(cosh(ϑ₀ᵢ₋ₖτᵢ))]
R∗ = 8.314510
R = 8.314472

τ = Tᵣ/T
δ = ρ/ρᵣ
(1/ρᵣ) = ∑∑xᵢxⱼβᵥ₋ᵢⱼγᵥ₋ᵢⱼ[(xᵢ+xⱼ)/(xᵢβᵥ₋ᵢⱼ^2 + xⱼ)]•1/8(1/∛ρcᵢ + 1/∛ρcⱼ)^2
Tᵣ = ∑∑xᵢxⱼβₜ₋ᵢⱼγₜ₋ᵢⱼ[(xᵢ+xⱼ)/(xᵢβₜ₋ᵢⱼ^2 + xⱼ)]•√(TcᵢTcⱼ)
aʳ = ∑xᵢaᵣᵢ(τ,δ) + ∑∑xᵢxⱼFᵢⱼaʳᵢⱼ(τ,δ)
aʳᵢ = ∑nᵢ₋ₖδ^(dᵢ₋ₖ)τ^(tᵢ₋ₖ)  + ∑nᵢ₋ₖδ^(dᵢ₋ₖ)τ^(tᵢ₋ₖ)exp(-δ^cᵢ₋ₖ)
aʳᵢⱼ = ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)τ^(tᵢⱼ₋ₖ)  + ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)τ^(tᵢⱼ₋ₖ)exp(ηᵢⱼ₋ₖ(δ-εᵢⱼ₋ₖ)^2 + βᵢⱼ₋ₖ(δ-γᵢⱼ₋ₖ)) 
```

## References

1. Kunz, O., & Wagner, W. (2012). The GERG-2008 wide-range equation of state for natural gases and other mixtures: An expansion of GERG-2004. Journal of Chemical and Engineering Data, 57(11), 3032–3091. [doi:10.1021/je300655b](https://doi.org/10.1021/je300655b)
"""
GERG2008