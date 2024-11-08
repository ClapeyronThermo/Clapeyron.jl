struct QuadraticDepartureParam <: EoSParam
    k0::PairParam{Float64}
    k1::PairParam{Float64}
end

@newmodelsimple QuadraticDeparture MultiFluidDepartureModel QuadraticDepartureParam

"""
    QuadraticDeparture <: MultiFluidDepartureModel
    QuadraticDeparture(components;
    userlocations = String[],
    verbose = false)

## Input parameters
- `k0`: Pair Parameter (`Float64`) - binary interaction parameter  (no units)
- `k1`: Pair Parameter (`Float64`) - binary, T-dependent interaction parameter `[K^-1]`

## Description

Departure that uses a quadratic mixing rule:

```
aᵣ = ∑xᵢxⱼaᵣᵢⱼ
aᵣᵢⱼ = 0.5*(aᵣᵢ + aᵣⱼ)*(1 - (k₀ + k₁T))
```

## References
1. Jäger, A., Breitkopf, C., & Richter, M. (2021). The representation of cross second virial coefficients by multifluid mixture models and other equations of state. Industrial & Engineering Chemistry Research, 60(25), 9286–9295. [doi:10.1021/acs.iecr.1c01186](https://doi.org/10.1021/acs.iecr.1c01186)
"""
QuadraticDeparture
default_locations(::Type{QuadraticDeparture}) = ["Empiric/departure/quadratic_departure_unlike.csv"]
default_references(::Type{QuadraticDeparture}) = ["10.1021/acs.iecr.1c01186","10.1016/j.fluid.2018.04.015"]
function transform_params(::Type{QuadraticDeparture},params,components)
    k0 = get(params,"k0",nothing)
    k1 = get(params,"k1",nothing)
    k0 === nothing && (k0 = PairParam("k0",components))
    k1 === nothing && (k1 = PairParam("k0",components))
    params["k0"] = k0
    params["k1"] = k1
    return params
end

function multiparameter_a_res(model,V,T,z,departure::QuadraticDeparture,δ,τ,∑z = sum(z))
    lnδ = log(δ)
    lnτ = log(τ)
    _0 = zero(lnδ+lnτ)
    n = length(model)
    aᵣₖ = fill(_0,length(model))
    m = model.pures
    Rinv = 1/Rgas(model)
    for i in 1:n
        mᵢ = m[i]
        aᵣᵢ[i] = reduced_a_res(mᵢ,δ,τ,lnδ,lnτ)*Rinv*Rgas(mᵢ)
    end
    k₀ = departure.params.k0
    k₁ = departure.params.k1
    aᵣ = _0
    for i in 1:n
        aᵢ = aᵣₖ[i]
        zᵢ = z[i]
        aᵣ += aᵢ*zᵢ*zᵢ
        for j in 1:(i-1)
            aᵢⱼ += 2*zᵢ*z[j]*0.5*(aᵢ + aᵣₖ[j])*(1 - k₀[i,j] - k₁[i,j]*T)
        end
    end
    return aᵣ/(∑z*∑z)
end

export QuadraticDeparture