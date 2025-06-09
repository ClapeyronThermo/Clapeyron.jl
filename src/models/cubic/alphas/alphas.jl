"""
    α_function(model::CubicModel,V,T,z,αmodel::AlphaModel)

Interface function used in cubic models. it should return a vector of αᵢ(T).

## Example:

```julia
function α_function(model::CubicModel,V,T,z,alpha_model::RKAlphaModel)
    return 1 ./ sqrt.(T ./ Tc)
end
```
"""
function α_function end

α_function(model,V,T,z) = α_function(model,V,T,z,model.alpha)
#all alphas at the moment don't have any need for recombine!
recombine_impl!(model::AlphaModel) = model
recombine_alpha!(model::CubicModel,alpha::AlphaModel) = recombine!(alpha)

struct SimpleAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

struct CPAAlphaParam <: EoSParam
    c1::SingleParam{Float64}
end

function fast_build_alpha(::Type{T}) where T <: AlphaModel
    if hasfield(T,:params)
        if fieldtype(T,:params) == SimpleAlphaParam
            return true
        elseif fieldtype(T,:params) == CPAAlphaParam
            return true
        end
    end
    return false
end

function __ignored_crit_params(alpha)
    if can_build_alpha_w(alpha)
        return ["Vc"]
    else
        return ["Vc","acentricfactor"]
    end
end

fast_build_alpha(T) = false

function init_alphamodel(alpha,components,params,userlocations = String[],verbose = [])
    #Base.Callable = Union{Type,Function}
    if alpha isa Base.Callable && alpha <: AlphaModel
        _comps = format_components(components)
        if fast_build_alpha(alpha) && isempty(userlocations)
            PARAM = fieldtype(alpha,:params)
            out_params = transform_params(PARAM,params,_comps)
            param = build_eosparam(PARAM,out_params)
            return alpha(_comps,param,default_references(typeof(alpha)))
        end
    end
    return init_model(alpha,components,userlocations,verbose)
end 

include("NoAlpha.jl")
include("ClausiusAlpha.jl")
include("RKAlpha.jl")
include("soave.jl")
include("soave2019.jl")
include("PRAlpha.jl")
include("Leibovici.jl")
include("PatelTejaAlpha.jl")
include("PTVAlpha.jl")
include("CPAAlpha.jl")
include("PR78Alpha.jl")
include("BM.jl")
include("Twu.jl")
include("MT.jl")
include("KUAlpha.jl")
include("RKPRAlpha.jl")

