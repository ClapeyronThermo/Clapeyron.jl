#=
abstract type ClausiusTranslationModel <: TranslationModel end

struct ClausiusTranslationParam <: EoSParam
    Vc::SingleParam{Float64}
end

@newmodelsimple ClausiusTranslation ClausiusTranslationModel ClausiusTranslationParam

"""
    ClausiusTranslation <: ClausiusTranslationModel

    ClausiusTranslation(components::Vector{String};
    userlocations::Vector{String}=String[],
    verbose::Bool=false)

## Input Parameters

- `vc`: Single Parameter (`Float64`) - Critical Volume `[m³/mol]`

## Description

Clausius Translation model for cubics:
```
V = V₀ + mixing_rule(cᵢ)
cᵢ = Zc_eosᵢ*R̄*Tci/Pcᵢ - Vcᵢ
Zc_eosᵢ = (1 - (Δ1+Δ2-1)*Ωbᵢ)/3
```
where `kᵢ` is the Cubic EoS calculated critical compresibility factor
## References

1. Clausius, R. (1880). Ueber das Verhalten der Kohlensäure in Bezug auf Druck, Volumen und Temperatur. Annalen der Physik, 245(3), 337–357. [doi:10.1002/andp.18802450302](https://doi.org/10.1002/andp.18802450302)

"""
ClausiusTranslation

export ClausiusTranslation
function ClausiusTranslation(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    Vc = params["Vc"]
    packagedparams = ClausiusTranslationParam(Vc)
    model = ClausiusTranslation(packagedparams, verbose=verbose)
    return model
end

function translation(model::CubicModel,V,T,z,translation_model::ClausiusTranslationModel)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = translation_model.params.Vc.values
    c = zeros(eltype(Tc),length(Tc))
    _,Ωb = ab_consts(model)
    Δ1,Δ2 = cubic_Δ(model,z)
    k = (1 - (Δ1+Δ2-1)*Ωb)/3
    if isone(length(Ωb))
        Ωbx = FillArrays.Fill(Ωb,length(z))
    else
        Ωbx = Ωb
    end        
    for i ∈ @comps
        Tci = Tc[i]
        Pci = Pc[i]
        ki = (1 + (Δ1+Δ2+1)*Ωbx[i])/3
        c[i] = ki*R̄*Tci/Pci - Vc[i]
    end
    return c
end=#