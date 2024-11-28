"""
    α_function(model::CubicModel,V,T,z,αmodel::AlphaModel)

Interface function used in cubic models. it should return a vector of αᵢ(T).

## Example:

```julia
function α_function!(α,model::CubicModel,alpha_model::RKAlphaModel)
    return 1 ./ sqrt.(T ./ Tc)
end
```
"""
function α_function end

α_function(model,V,T,z) = α_function(model,V,T,z,model.alpha)

#to allow dispatch for single comp in old versions
α_function(model,V,T) = α_function(model,V,T,SA[1.0],model.alpha)

#In the future, we could use the optional allocation argument to stop allocating alphas
α_function(model,V,T,z,alphamodel) = _α_function(model,V,T,z,alphamodel)

function α_function(model,V,T,z,alphamodel)
    _0 = zero(Base.promote_eltype(model,T))
    nc = length(model)
    if nc <= 8
        #view of a static array of constant length
        #this allows non-allocating alpha functions.
        αtup = ntuple(i -> i <= nc ? α_function(model,alphamodel,T,i) : _0,Val{8}())
        return view(SVector(αtup),1:nc)
    else
        αvec = fill(_0,nc)
        α_function!(αvec,model,alphamodel,T)
        return αvec
    end
end

#all alphas at the moment don't have any need for recombine!
recombine_impl!(model::AlphaModel) = model
recombine_alpha!(model::CubicModel,alpha::AlphaModel) = recombine!(alpha)

struct SimpleAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

function can_build_alpha_w(::Type{T}) where T <: AlphaModel
    if hasfield(T,:params)
        if fieldtype(T,:params) == SimpleAlphaParam
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

can_build_alpha_w(T) = false

function init_alphamodel(alpha,components,w = nothing,userlocations = String[],verbose = [])
    #Base.Callable = Union{Type,Function}
    if alpha isa Base.Callable && alpha <: AlphaModel
        if can_build_alpha_w(alpha) && w !== nothing && isempty(userlocations)
            param = SimpleAlphaParam(w)
            return alpha(format_components(components),param,default_references(typeof(alpha)))
        end
    end
    return init_model(alpha,components,userlocations,verbose)
end

function α_function!(vec,model,alphamodel,T)
    @inbounds for i in 1:length(model)
        vec[i] = α_function(model,alphamodel,T,i)
    end
    return vec
end

include("NoAlpha.jl")
include("ClausiusAlpha.jl")
include("RKAlpha.jl")
include("soave.jl")
include("soave2019.jl")
include("PRAlpha.jl")
include("PatelTejaAlpha.jl")
include("PTVAlpha.jl")
include("CPAAlpha.jl")
include("sCPAAlpha.jl")
include("PR78Alpha.jl")
include("BM.jl")
include("Twu.jl")
include("MT.jl")
include("KUAlpha.jl")
include("RKPRAlpha.jl")

