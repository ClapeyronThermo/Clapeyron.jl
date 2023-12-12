struct SLKRuleParam <: EoSParam
    k::PairParam{Float64}
end

@newmodelsimple SLKRule SLMixingRule SLKRuleParam

default_locations(::Type{SLKRule}) = ["LatticeFluid/SanchezLacombe/mixing/k0k1l_unlike.csv"]
default_ignore_missing_singleparams(::Type{SLKRule}) = ["k"]

function transform_params(::Type{SLKRule},params,components)
    n = length(components)

    if haskey(params,"k0") && !haskey(params,"k")
        params["k"] = params["k0"]
    end
    
    k0 = get!(params,"k") do
        PairParam("k",components,zeros(n,n))
    end

    return params
end

export SLKRule

"""
    SLKRule(components; userlocations=String[], verbose=false)
     
## Input parameters
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)

Constant Kᵢⱼ mixing rule for Sanchez-Lacombe:

```
εᵢⱼ = √εᵢεⱼ*(1-kᵢⱼ)
vᵢⱼ = (vᵢ + vⱼ)/2
ϕᵢ = rᵢ*xᵢ/r̄
εᵣ = ΣΣϕᵢϕⱼεᵢⱼ
vᵣ = ΣΣϕᵢϕⱼvᵢⱼ
```

"""
SLKRule

function set_k!(model::SanchezLacombe{SLKRule},k)
    check_arraysize(model,k)
    model.mixing.k.values = k
end

function __SL_get_k(model::Clapeyron.SanchezLacombe,mixing::Clapeyron.SLKRule)
    return copy(mixing.k.values)
end


function sl_mix(unmixed_vol,unmixed_epsilon,mixmodel::SLKRule)
    #dont mind the function names, it performs the correct mixing
    premixed_vol = epsilon_LorentzBerthelot(unmixed_vol)
    premixed_epsilon = sigma_LorentzBerthelot(unmixed_epsilon)
    return premixed_vol,premixed_epsilon
end

function sl_mix!(unmixed_vol,unmixed_epsilon,mixmodel::SLKRule)
    #dont mind the function names, it performs the correct mixing
    epsilon_LorentzBerthelot!(unmixed_vol)
    sigma_LorentzBerthelot!(unmixed_epsilon)
end

function mix_vε(model::SanchezLacombe,V,T,z,mix::SLKRule,r̄,Σz = sum(z))
    v = model.params.vol.values
    ε = model.params.epsilon.values
    isone(length(z)) && return (only(v),only(ε))
    r =  model.params.segment.values
    k = mix.params.k.values
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
            εij = ε[i,j]*(1-k[i,j])
            ε_r += ϕiϕj*εij
        end
    end
    return v_r,ε_r
end
