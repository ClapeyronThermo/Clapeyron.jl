struct SLk0k1lMixingRule <: SLMixingRule
    components::Vector{String}
    k0::PairParam{Float64}
    k1::PairParam{Float64}
    l::PairParam{Float64}
end

@registermodel SLk0k1lMixingRule

function SLk0k1lMixingRule(components; userlocations=String[], verbose=false)
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
    r =  model.params.segment.values
    v = model.params.vol.values
    k0 = mix.k0.values
    k1 = mix.k1.values
    ε = model.params.epsilon.values
    r̄inv = one(r̄)/r̄
     ϕ = @. r* z* r̄inv
    v_r = zero(V+T+first(z))
    ε_r = v_r
    Σz2 = 1/(Σz*Σz)
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
    return v_r*Σz2,ε_r*Σz2
end