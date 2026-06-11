module StaticForwardDiffTag

using ForwardDiff
using ForwardDiff: Tag, Dual, tagcount
import ForwardDiff: ≺

"""
    StaticTag{F,V}

A replacement of `ForwardDiff.Tag{F,V}` that supports a faster nested comparison.
"""
struct StaticTag{F,V} end

function StaticTag(f::F, ::Type{V}) where {F,V}
    tagcount(Tag{F,V}) # trigger generated function, but with the tag type.
    StaticTag{F,V}()
end

function StaticTag(f::Type{F}, ::Type{V}) where {F,V}
    tagcount(Tag{F,V}) # trigger generated function, but with the tag type.
    StaticTag{F,V}()
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

@inline function ≺(::Type{StaticTag{F1,V1}}, ::Type{Tag{F2,V2}}) where {F1,V1,F2,V2}
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

@inline function ≺(::Type{Tag{F1,V1}}, ::Type{StaticTag{F2,V2}}) where {F1,V1,F2,V2}
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

@inline function ≺(::Type{StaticTag{F1,V1}}, ::Type{StaticTag{F2,V2}}) where {F1,V1,F2,V2}
    T1 = Tag{F1,V1}
    T2 = Tag{F2,V2}
    d1 = tagdepth(T1)
    d2 = tagdepth(T2)
    if d1 != d2
        genuinely_nested = d1 > d2 ? contains_tag(T1, T2) : contains_tag(T2, T1)
        genuinely_nested && return d1 < d2   # depth wins ONLY if nesting is real
    end
    throw(DualMismatchError(T1,T2))
end


ForwardDiff.checktag(::Type{StaticTag{FT,VT}}, f::F, x::AbstractArray{V}) where {FT,VT,F,V} =
    throw(ForwardDiff.InvalidTagException{Tag{F,V},Tag{FT,VT}}())

ForwardDiff.checktag(::Type{StaticTag{F,V}}, f::F, x::AbstractArray{V}) where {F,V} = true

# no easy way to check Jacobian tag used with Hessians as multiple functions may be used
ForwardDiff.checktag(::Type{StaticTag{FT,VT}}, f::F, x::AbstractArray{V}) where {FT<:Tuple,VT,F,V} = true

"""
    ∂Deferred{F,P}

A deferred evaluation system to skip the Closure perturbation confusion problem in ForwardDiff.jl
"""
struct ∂Deferred{T,V,F,P}
    f::F
    p::P
end

"""
    ∂Tag{T}

a static tag type, used with deferred evaluation
"""
struct ∂Tag{T} end


∂Deferred(f::F,p::P) where {F,P} = ∂Deferred(f,p,∂Tag{F})

∂Deferred(f::F,p::V,::T) where {F,V <: Number,T} = ∂Deferred{T,V,F,V}(f,p)
∂Deferred(f::F,p::AbstractArray{V},::T) where {F,V <: Number,T} = ∂Deferred{T,V,F,typeof(p)}(f,p)

function ∂Deferred(f::F,p::Tup,::T) where {F,Tup<:Tuple,T}
    V = Base.promote_eltype(p...)
    return ∂Deferred{T,V,F,Tup}(f,p)
end

(∂::∂Deferred{F,P})(x) where {F,P} = ∂.f(∂.p)(x)

function ForwardDiff.Tag(f::∂Deferred{T,V1,F,P},::Type{V2}) where {T,V1,F,P,V2}
    return StaticTag(T,deferred_valtype(V1,V2))
end

function StaticTag(f::∂Deferred{T,V1,F,P},::Type{V2}) where {T,V1,F,P,V2}
    return StaticTag(T,deferred_valtype(V1,V2))
end

function ForwardDiff.Tag(::Type{∂Deferred{T,V1,F,P}},::Type{V2}) where {T,V1,F,P,V2}
    return StaticTag(T,deferred_valtype(V1,V2))
end

function StaticTag(::Type{∂Deferred{T,V1,F,P}},::Type{V2}) where {T,V1,F,P,V2}
    return StaticTag(T,deferred_valtype(V1,V2))
end


deferred_valtype(::Type{V},::Type{T}) where {V,T} = promote_type(V,T)

end #module

