struct SLKRule <: SLMixingRule
    components::Vector{String}
    k::PairParam{Float64}
end

@registermodel SLKRule

function sl_mix(unmixed_vol,unmixed_epsilon,mixmodel::SLKRule)
    #dont mind the function names, it performs the correct mixing
    premixed_vol= epsilon_LorentzBerthelot(unmixed_vol)
    premixed_epsilon = sigma_LorentzBerthelot(unmixed_epsilon)
    return premixed_vol,premixed_epsilon
end

function SLKRule(components; userlocations=String[], verbose=false)
    params = getparams(components, ["LatticeFluid/SanchezLacombe/mixing/k0k1l_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k = params["k0"]
    model = SLKRule(components,k)
    return model
end

function mix_vε(model::SanchezLacombe,V,T,z,mix::SLKRule,r̄,Σz)
    v = model.params.vol.values
    ε = model.params.epsilon.values
    isone(length(z)) && return (only(v),only(ε))
    r =  model.params.segment.values
    k = mix.k.values
    r̄inv = one(r̄)/r̄
    ϕ = @. r* z* r̄inv/Σz
    v_r = zero(V+T+first(z))
    ε_r = v_r
    Σz2 = 1/(Σz*Σz)
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