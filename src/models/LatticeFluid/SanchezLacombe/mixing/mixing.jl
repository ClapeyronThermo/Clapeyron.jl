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

function init_slmixing(model::EoSModel,components,userlocations,mixing_userlocations,verbose)
    return model
end

function init_slmixing(model,components,params,mixing_userlocations,verbose)
    if any(z -> haskey(params,z),("k","k0","k1","l")) && model <: SLMixingRule
        paramstype = fieldtype(model,:params)
        pnew = transform_params(model,params,components)
        mixing_params = paramstype(pnew)
        if verbose
            @info "Building an instance of $(info_color(string(model))) with components $components"
        end
        return model(components,mixing_params,default_references(model))
    else
        return init_model(model,components,mixing_userlocations,verbose)
    end
end