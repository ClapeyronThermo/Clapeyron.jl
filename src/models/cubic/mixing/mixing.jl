#this is used for some dispatches
abstract type ActivityMixingRule <: MixingRule end
"""
    mixing_rule(model::CubicModel,V,T,z,mixing_model::MixingRule,α,a,b,c)

Interface function used by cubic models. with matrices `a` and `b`, vectors `α` and `c`, a `model::CubicModel` and `mixing_model::MixingRule`, returns
the scalars `ā`,`b̄` and `c̄`, corresponding to the values mixed by the amount of components and the specifics of the mixing rule.

## Example
```julia
function mixing_rule(model::CubicModel,V,T,z,mixing_model::vdW1fRule,α,a,b,c)
    ∑z = sum(z)
    ā = dot(z .* sqrt(α),a,z .* sqrt(α))/(∑z*∑z) #∑∑aᵢⱼxᵢxⱼ√(αᵢαⱼ)
    b̄ = dot(z,b,z)/(∑z*∑z)  #∑∑bᵢⱼxᵢxⱼ
    c̄ = dot(z,c)/∑z ∑cᵢxᵢ
    return ā,b̄,c̄
end
```
"""
function mixing_rule end

function init_model(model::MixingRule,components,activity,userlocations,activity_userlocations,verbose)
    return model
end

function init_model(model::Type{<:MixingRule},components,activity,userlocations,activity_userlocations,verbose)
    if verbose
        @info "Building an instance of $(info_color(string(model))) with components $components"
    end
    return model(components;activity,userlocations,activity_userlocations,verbose)
end

#used in CompositeModel.jl
init_mixing_act = init_model_act

function infinite_pressure_gibbs_correction(model::DeltaCubicModel,T,z)
    Δ1,Δ2 = cubic_ΔT(model,T,z)
    if Δ1==Δ2
        return 1/(1-Δ1)
    else
        return -log((1-Δ1)/(1-Δ2))/(Δ1 - Δ2)
    end
end

function infinite_pressure_gibbs_correction(model::vdWModel,T,z)
    return -1.0
end

#default
function recombine_mixing!(model,mixing_model,k = nothing, l = nothing)
    recombine!(mixing_model)
    c = c_premixing(model)
    ab_premixing(model,mixing_model,k,l)
    return mixing_model
end

function recombine_impl!(model::ActivityMixingRule)
    recombine!(model.activity)
    return model
end



include("vdW1f.jl")
include("Kay.jl")
include("HV.jl")
include("MHV1.jl")
include("MHV2.jl")
include("LCVM.jl")
include("WS.jl")
include("modWS.jl")
include("PSRK.jl")
include("VTPR.jl")
include("UMR.jl")
include("QCPR.jl")
include("PPR78.jl")
include("gEr.jl")