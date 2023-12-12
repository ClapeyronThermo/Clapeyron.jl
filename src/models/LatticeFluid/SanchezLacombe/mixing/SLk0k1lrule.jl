struct SLk0k1lMixingRuleParam <: EoSParam
    k0::PairParam{Float64}
    k1::PairParam{Float64}
    l::PairParam{Float64}
end

@newmodelsimple SLk0k1lMixingRule SLMixingRule SLk0k1lMixingRuleParam


default_locations(::Type{SLk0k1lMixingRule}) = ["LatticeFluid/SanchezLacombe/mixing/k0k1l_unlike.csv"]
default_ignore_missing_singleparams(::Type{SLk0k1lMixingRule}) = ["k","k0","k1","l"]

function transform_params(::Type{SLk0k1lMixingRule},params,components)
    if haskey(params,"k") && !haskey(params,"k0")
        params["k0"] = params["k"]
    end

    n = length(components)

    k0 = get!(params,"k0") do
        PairParam("k0",components,zeros(n,n))
    end

    k1 = get!(params,"k1") do
        PairParam("k1",components,zeros(n,n))
    end

    l = get!(params,"l") do
        PairParam("l",components,zeros(n,n))
    end
    return params
end

export SLk0k1lMixingRule

"""
    SLKRule(components; userlocations=String[], verbose=false)

## Input parameters
- `k0`,`k`: Pair Parameter (`Float64`, optional) - Binary Interaction Parameter (no units)
- `k1`: Pair Parameter (`Float64`, optional) - Binary Interaction Parameter (no units)
- `l`: Pair Parameter (`Float64`,optional) - Binary Interaction Parameter (no units)

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

function sl_mix(unmixed_vol,unmixed_epsilon,mixmodel::SLk0k1lMixingRule)
    #dont mind the function names, it performs the correct mixing
    premixed_vol = epsilon_LorentzBerthelot(unmixed_vol,mixmodel.l)
    premixed_epsilon = sigma_LorentzBerthelot(unmixed_epsilon)
    return premixed_vol,premixed_epsilon
end

function sl_mix!(unmixed_vol,unmixed_epsilon,mixmodel::SLk0k1lMixingRule)
    epsilon_LorentzBerthelot!(unmixed_vol,mixmodel.l)
    sigma_LorentzBerthelot!(unmixed_epsilon)
end

function __SL_get_k(model::SanchezLacombe,mixing::SLk0k1lMixingRule)
    return copy(mixing.k0.values),copy(mixing.k1.values)
end

function set_k!(model::SanchezLacombe{SLk0k1lMixingRule},k0,k1)
    check_arraysize(model,k0)
    check_arraysize(model,k1)
    model.mixing.k0.values = k0
    model.mixing.k1.values = k1
end

function set_k!(model::SanchezLacombe{SLk0k1lMixingRule},k0)
    check_arraysize(model,k0)
    model.mixing.k0.values .= k0
    n = length(model)
    model.mixing.k1.values .= FillArrays.Zeros(n,n)
end

function __SL_get_k(model::SanchezLacombe,mixing::SLKRule)
    return copy(mixing.k.values)
end

function __SL_get_l(model::SanchezLacombe,mixing::SLk0k1lMixingRule)
    return copy(mixing.l.values) 
end

function set_l!(model::SanchezLacombe{SLk0k1lMixingRule},l)
    l0 = model.mixing.l.values
    epsilon_LorentzBerthelot!(model.params.vol,l)
    l0 .= l
end

function mix_vε(model::SanchezLacombe,V,T,z,mix::SLk0k1lMixingRule,r̄,Σz)
    ε = model.params.epsilon.values
    v = model.params.vol.values
    isone(length(z)) && return (only(v),only(ε))
    r =  model.params.segment.values
    k0 = mix.params.k0.values
    k1 = mix.params.k1.values
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
            kij = k0[i,j] + (1 - δij) * (dot(ϕ,kmi) + dot(ϕ,kmj))
            εij = ε[i,j]*(1 - kij)
            ε_r += ϕiϕj*εij
        end
    end
    return v_r,ε_r
end
