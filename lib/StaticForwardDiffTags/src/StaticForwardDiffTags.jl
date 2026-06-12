module StaticForwardDiffTags

using ForwardDiff
using ForwardDiff: Tag, Dual, tagcount, checktag, Chunk
using ForwardDiff.DiffResults
import ForwardDiff: ≺

"""
    STag{F,V}

A replacement of `ForwardDiff.Tag{F,V}` that supports a faster nested comparison.
"""
struct STag{F,V} end

"""
    ∂Tag{F}

a static "function type",marking that the function parametrized is in fact pure, so no `ForwardDiff.tagcount` is needed at evaluation time
"""
struct ∂Tag{F} end

@inline function STag(f::F, ::Type{V}) where {F,V}
    tagcount(Tag{F,V}) # trigger generated function, but with the tag type.
    STag{F,V}()
end

@inline function STag(f::Type{F}, ::Type{V}) where {F,V}
    tagcount(Tag{F,V}) # trigger generated function, but with the tag type.
    STag{F,V}()
end

@inline tagdepth(::Type) = 0
@inline tagdepth(::Type{<:Dual{T,V,N}}) where {T,V,N} = 1 + tagdepth(V)
@inline tagdepth(::Type{<:Tag{F,V}}) where {F,V} = 1 + tagdepth(V)

# Does Tag `target` genuinely appear in X's value-type / tag chain?
# Deliberately does NOT recurse closure type params F.
@inline contains_tag(::Type, target) = false
@inline contains_tag(::Type{<:Tag{F,V}}, target) where {F,V} = contains_tag(V, target)
@inline contains_tag(::Type{<:Dual{T,V,N}}, target) where {T,V,N} =
    (T === target) || contains_tag(T, target) || contains_tag(V, target)

#dynamic fallbacks, with jarod ForwardDiff#807 patch
@inline function ≺(::Type{STag{F1,V1}}, ::Type{Tag{F2,V2}}) where {F1,V1,F2,V2}
    T1 = Tag{F1,V1}
    T2 = Tag{F2,V2}
    d1 = tagdepth(T1)
    d2 = tagdepth(T2)
    if d1 != d2
        genuinely_nested = d1 > d2 ? contains_tag(T1, T2) : contains_tag(T2, T1)
        genuinely_nested && return d1 < d2   # depth wins ONLY if nesting is real
    end
    return tagcount(T1) < tagcount(T2)
end

@inline function ≺(::Type{Tag{F1,V1}}, ::Type{STag{F2,V2}}) where {F1,V1,F2,V2}
    T1 = Tag{F1,V1}
    T2 = Tag{F2,V2}
    d1 = tagdepth(T1)
    d2 = tagdepth(T2)
    if d1 != d2
        genuinely_nested = d1 > d2 ? contains_tag(T1, T2) : contains_tag(T2, T1)
        genuinely_nested && return d1 < d2   # depth wins ONLY if nesting is real
    end
    return tagcount(T1) < tagcount(T2)
end

@inline is_pure_f(::Type{F}) where {F} = false
@inline is_pure_f(::Type{F}) where {F <: ∂Tag} = true

#Static checking of STags
@inline function ≺(::Type{STag{F1,V1}}, ::Type{STag{F2,V2}}) where {F1,V1,F2,V2}
    T1 = Tag{F1,V1}
    T2 = Tag{F2,V2}
    d1 = tagdepth(T1)
    d2 = tagdepth(T2)
    if d1 != d2
        #only deferred methods being used: depth is direct indicator of nesting
        is_pure_f(F1) && is_pure_f(F2) && return d1 < d2
        genuinely_nested = d1 > d2 ? contains_tag(T1, T2) : contains_tag(T2, T1)
        genuinely_nested && return d1 < d2   # depth wins ONLY if nesting is real
    end
    throw(ForwardDiff.DualMismatchError(T1,T2))
end

ForwardDiff.checktag(::Type{STag{FT,VT}}, f::F, x::AbstractArray{V}) where {FT,VT,F,V} =
    throw(ForwardDiff.InvalidTagException{Tag{F,V},Tag{FT,VT}}())

ForwardDiff.checktag(::Type{STag{F,V}}, f::F, x::AbstractArray{V}) where {F,V} = true

# no easy way to check Jacobian tag used with Hessians as multiple functions may be used
ForwardDiff.checktag(::Type{STag{FT,VT}}, f::F, x::AbstractArray{V}) where {FT<:Tuple,VT,F,V} = true

"""
    SDiffFunction{F,P}

A deferred evaluation system to skip the Closure perturbation confusion problem in ForwardDiff.jl.
"""
struct SDiffFunction{T,V,F,P}
    f::F
    p::P
end

SDiffFunction(f::F,p::P) where {F,P} = SDiffFunction(f,p,∂Tag{F})

SDiffFunction(f::F,p::V,::T) where {F,V <: Number,T} = SDiffFunction{T,V,F,V}(f,p)
SDiffFunction(f::F,p::AbstractArray{V},::T) where {F,V <: Number,T} = SDiffFunction{T,V,F,typeof(p)}(f,p)

function SDiffFunction(f::F,p::Tup,::T) where {F,Tup<:Tuple,T}
    V = Base.promote_eltype(p...)
    return SDiffFunction{T,V,F,Tup}(f,p)
end

(∂::SDiffFunction{F,P})(x) where {F,P} = ∂.f(∂.p)(x)
(∂::SDiffFunction{F,P})(x,y) where {F,P} = ∂.f(∂.p)(x,y)

#pure function api
(∂::SDiffFunction{F,Nothing})(x,y) where {F} = ∂.f(x,y)
(∂::SDiffFunction{F,Nothing})(x) where {F} = ∂.f(x)

@inline function STag(f::SDiffFunction{T,V1,F,P},::Type{V2}) where {T,V1,F,P,V2}
    return STag(T,deferred_valtype(V1,V2))
end

@inline function STag(::Type{SDiffFunction{T,V1,F,P}},::Type{V2}) where {T,V1,F,P,V2}
    return STag(T,deferred_valtype(V1,V2))
end

@inline deferred_valtype(::Type{V1},::Type{V2}) where {V1,V2} = promote_type(V1,V2)

"""
    maketag(f,v)
    maketag(f,::Type{V})

returns a ForwardDiff-Compatible tag type. For generic functions it dispatches to `ForwardDiff.Tag`.
"""
maketag(f::F,v::V) where {F,V} = Tag(f,V)
maketag(f::F,::Type{V}) where {F,V} = Tag(f,V)
maketag(f::F,v::V) where {F<:SDiffFunction,V} = STag(f,V)
maketag(f::F,v::Type{V}) where {F<:SDiffFunction,V} = STag(f,V)

@inline maketagtype(f::F,v::V) where {F,V} = typeof(maketag(f,v))

include("config.jl") #config overloads
include("derivative.jl") #derivative overloads

end #module

