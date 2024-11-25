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

function _α_function(model,V,T,z,alphamodel)
    return LazyAlphaVector(model,alphamodel,T)
    α = zeros(Base.promote_eltype(model,T),length(model))
    α_function!(α,model,alphamodel,T)
    return α
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

mutable struct LazyAlphaVector{M,α,TT} <: AbstractVector{TT}
    @const model::M
    @const alphamodel::α
    @const T::TT
    @const allocated::Bool
    @const nc::Int
    vector::Vector{TT}
    function LazyAlphaVector(model::M,alpha::α,T::T1) where {M,α,T1}
        TT = Base.promote_eltype(model,T)
        nmax = alpha_lazy_limit(model)
        nc = length(model)
        needs_alloc = nc > nmax
        αvec = new{M,α,TT}(model,alpha,convert(TT,T),needs_alloc,nc)
        if needs_alloc
            vec = Vector{TT}(undef,nc)
            α_function!(vec,model,alpha,T)
            αvec.vector = vec
        end
        return αvec
    end
end

alpha_lazy_limit(model::CubicModel) = min(alpha_lazy_limit(model.mixing),alpha_lazy_limit(model.alpha))
alpha_lazy_limit(model) = 5

function α_function!(vec,model,alphamodel,T)
    for i in 1:length(model)
        vec[i] = α_function(model,alphamodel,T,i)
    end
    return vec
end

Base.size(x::LazyAlphaVector) = (x.nc,)
Base.length(x::LazyAlphaVector) = x.nc
Base.eltype(x::LazyAlphaVector{M,α,TT}) where {M,α,TT} = TT
Base.@propagate_inbounds function Base.getindex(x::LazyAlphaVector{M,α,TT},i::Int)::TT where {M,α,TT}
    if x.allocated
        @inbounds return x.vector[i]::TT
    else
        return α_function(x.model,x.alphamodel,x.T,i)::TT
    end
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

