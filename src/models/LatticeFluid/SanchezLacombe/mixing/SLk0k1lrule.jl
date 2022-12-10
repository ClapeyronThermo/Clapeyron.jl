struct SLk0k1lMixingRule <: SLMixingRule
    components::Vector{String}
    k0::PairParam{Float64}
    k1::PairParam{Float64}
    l::PairParam{Float64}
end

@registermodel SLk0k1lMixingRule
export SLk0k1lMixingRule

"""
    SLKRule(components; userlocations=String[], verbose=false)
     
## Input parameters
- `k0`: Pair Parameter (`Float64`) - Binary Interaction Parameter (no units)
- `k1`: Pair Parameter (`Float64`) - Binary Interaction Parameter (no units)
- `l`: Pair Parameter (`Float64`) - Binary Interaction Parameter (no units)

Neau's Consistent k₀,k₁,l mixing rule for Sanchez-Lacombe:

```
εᵢⱼ = √εᵢεⱼ
vᵢⱼ = (1 - lᵢⱼ)(vᵢ + vⱼ)/2
ϕᵢ = rᵢ*xᵢ/r̄
εᵣ = ΣΣϕᵢϕⱼεᵢⱼ*(1 - k₀ᵢⱼ + (1 - δᵢⱼ)(Σϕₖk₁ᵢₖ + Σϕₖk₁ₖⱼ))
vᵣ = ΣΣϕᵢϕⱼvᵢⱼ
```
Where `δᵢⱼ` is `i == j ? 1 : 0`

## References
1. Neau, E. (2002). A consistent method for phase equilibrium calculation using the Sanchez–Lacombe lattice–fluid equation-of-state. Fluid Phase Equilibria, 203(1–2), 133–140. [doi:10.1016/s0378-3812(02)00176-0](https://doi.org/10.1016/s0378-3812(02)00176-0)
"""
SLk0k1lMixingRule

function SLk0k1lMixingRule(components; userlocations=String[], verbose=false, kwargs...)
    params = getparams(components, ["LatticeFluid/SanchezLacombe/mixing/k0k1l_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k0 = params["k0"]
    k1 = params["k1"]
    l = params["l"]
   model = SLk0k1lMixingRule(components,k0,k1,l)
    return model
end

function sl_mix(unmixed_vol,unmixed_epsilon,mixmodel::SLk0k1lMixingRule)
    #dont mind the function names, it performs the correct mixing
    premixed_vol= epsilon_LorentzBerthelot(unmixed_vol,mixmodel.l)
    premixed_epsilon = sigma_LorentzBerthelot(unmixed_epsilon)
    return premixed_vol,premixed_epsilon
end

function mix_vε(model::SanchezLacombe,V,T,z,mix::SLk0k1lMixingRule,r̄,Σz)
    ε = model.params.epsilon.values
    v = model.params.vol.values
    isone(length(z)) && return (only(v),only(ε))
    r =  model.params.segment.values
    k0 = mix.k0.values
    k1 = mix.k1.values
    r̄inv = one(r̄)/r̄
    ϕ = @. r* z* r̄inv/Σz
    v_r = zero(V+T+first(z))
    ε_r = v_r
    for i in @comps
        for j in @comps
            ϕi = ϕ[i]
            ϕj = ϕ[j]
            ϕiϕj = ϕi*ϕj 
            v_r += ϕiϕj*v[i,j]
            δij = (i===j)
            kmi = view(k1,:,i)
            kmj = view(k1,:,j)
            kij = k0[i,j] + (1-δij)* (dot(ϕ,kmi) + dot(ϕ,kmj))
      
            εij = ε[i,j]*(1-kij)
            ε_r += ϕiϕj*εij
        end
    end
    return v_r,ε_r
end
