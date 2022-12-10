abstract type SLMixingRule <: EoSModel end


"""
    mix_vε(model::SanchezLacombeModel,V,T,z,mix::SLMixingRule,r̄ = @f(rmix),∑z = sum(z))

Function used to dispatch on the different mixing rules available for Sanchez-Lacombe.

## Example:
```julia
function mix_vε(model::SanchezLacombe,V,T,z,mix::SLKRule,r̄,Σz = sum(z))
    v = model.params.vol.values
    ε = model.params.epsilon.values
    r =  model.params.segment.values
    k = mix.k.values
    x = z ./ Σz
    ϕ = @. r * x / r̄
    εᵣ = sum(ε[i,j]*(1-k[i,j])*ϕ[i]*ϕ[j] for i ∈ @comps for j ∈ @comps)
    vᵣ = sum(v[i,j]*ϕ[i]*ϕ[j] for i ∈ @comps for j ∈ @comps)
    return vᵣ,εᵣ
```
"""
function mix_vε end