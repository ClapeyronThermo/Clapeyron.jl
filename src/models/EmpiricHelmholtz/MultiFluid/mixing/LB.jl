"""
    LorentzBerthelotMixing::AsymmetricMixing
    LorentzBerthelotMixing(components;
    userlocations = String[],
    verbose = false)

## Input parameters
- `k`: Pair Parameter (`Float64`) - binary interaction parameter for temperature  (no units)
- `l`: Pair Parameter (`Float64`) - binary interaction parameter for volume (no units)

## Description
Lorentz-Berthelot Mixing for MultiParameter EoS models:

```
τ = T̄/T
δ = V̄/V
V̄ = ∑xᵢxⱼ * Vᵣᵢⱼ * (1 - lᵢⱼ)
T̄ = ∑xᵢxⱼ * Tᵣᵢⱼ * (1 - kᵢⱼ)
Vᵣᵢⱼ = 0.125*(∛Vcᵢ + ∛Vcⱼ)^3
Tᵣᵢⱼ = √(Tcᵢ*Tcⱼ)
```
missing parameters will be assumed `kᵢⱼ = lᵢⱼ = 0`
"""
function LorentzBerthelotMixing(components;userlocations = String[],verbose = false)
    params = getparams(components,["Empiric/mixing/LorentzBerthelotMixing/LorentzBerthelotMixing_unlike.csv"]; userlocations = userlocations, verbose = verbose)
    l = get(params,"l",nothing)
    k = get(params,"k",nothing)
    k === nothing && (k = PairParam("gamma_T",components))
    l === nothing && (l = PairParam("gamma_v",components))
    k .= 1 .- k #invert
    l .= 1 .- l
    beta_v = PairParam("beta_v",components)
    beta_T = PairParam("beta_T",components)
    beta_v .= 1
    beta_T .= 1
    pkgparams = AsymmetricMixingParam(k,l,beta_T,beta_v)
    n = length(components)
    k.ismissingvalues .= false
    l.ismissingvalues .= false
    beta_v.ismissingvalues .= false
    beta_T.ismissingvalues .= false
    references = ["Klimeck, Ph.D. dissertation"]
    return AsymmetricMixing(pkgparams,verbose = verbose,references = references)
end

export LorentzBerthelotMixing