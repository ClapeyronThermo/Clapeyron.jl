abstract type MixingRule <:EoSModel end



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
    verbose && @info("""Now creating mixing model:
    $model""")
    return model(components;activity,userlocations,activity_userlocations,verbose)
end

include("vdW1f.jl")
include("Kay.jl")
include("HV.jl")
include("MHV1.jl")
include("MHV2.jl")
include("LCVM.jl")
include("WS.jl")
include("PSRK.jl")
include("VTPR.jl")
include("UMR.jl")