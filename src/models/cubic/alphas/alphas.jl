"""
    α_function(model::CubicModel,V,T,z,αmodel::AlphaModel)

Interface function used in cubic models. It should return a vector of αᵢ(T).

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

function fast_build_alpha(T,params)
    if hasfield(T,:params)
        P = fieldtype(T,:params)
        if P == SimpleAlphaParam
            haskey(params,"acentricfactor") || return false
            w = params["acentricfactor"]
            return !any(w.ismissingvalues)
        elseif P == CPAAlphaParam
            haskey(params,"c1") || return false
            c1 = params["c1"]
            return !any(c1.ismissingvalues)
        end
    end
    return false
end

function init_alphamodel(alpha,components,params,userlocations = String[],verbose = [])
    #Base.Callable = Union{Type,Function}
    if alpha isa Base.Callable && alpha <: AlphaModel
        _comps = format_components(components)
        if fast_build_alpha(alpha,params) && isempty(userlocations)
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
include("MathiasCopemanAlpha.jl")
