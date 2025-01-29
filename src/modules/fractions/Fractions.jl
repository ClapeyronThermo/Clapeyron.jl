#= Composition inverses and mixes
a composition can be seen as a point in a Simplex manifold.
for that reason, you can define some common operations that are useful.
the
=#
module Fractions

    import ..Solvers
    frac(x) = x ./ Base.sum(x)
    sum(x1,x2) = frac(x1 .* x2)
    mul(x1,λ)  = frac(x1 .^λ)
    neg(x1) = frac(1 ./ x1)
    zeros(T,n) = fill(one(T)/n,n)
    zeros(n) = zeros(Float64,n)

    struct FractionVector{T,V} <: AbstractVector{T}
        vec::V ##Vector storing j != i values
        val::T #value at i
        i::Int #index i
    end
    #=

    Fraction Vector
    useful when expressing fractions in n-1 components.
    the last component is calculated at build time.
    it allocates less than creating a new vector or appending.
    =#
    ##
    function FractionVector(v::AbstractVector,i = length(v) + 1)
        a = Base.one(eltype(v))
        # any(x->x<0,v) && throw(DomainError(v,"all elements of a fraction vector should be positive."))
        a -=Base.sum(v)
        # a < 0 && throw(DomainError(a,"the values of the input vector add to more than one"))
        return FractionVector(v,a,i)
    end

    function FractionVector(v::Real,i = 2)
        a = Base.one(v)
        # (v < zero(v)) && throw(DomainError(v,"all elements of a fraction vector should be positive."))
        a -= v
        # a < 0 && throw(DomainError(a,"the values of the input vector add to more than one"))
        return FractionVector(v,a,i)
    end


    @inline Base.eltype(v::FractionVector{T}) where T = T
    @inline Base.length(v::FractionVector)::Int = Base.length(v.vec) + 1

    @inline function Base.length(v::FractionVector{T,<:Real})::Int where T
        return 2
    end

    @inline Base.size(v::FractionVector) = (length(v),)
    @inline function Base.getindex(v::FractionVector,j::Int)
        @boundscheck checkbounds(v, j)
        i = v.i
        j == i && return v.val
        ii = ifelse(j < i,j,j-1)
        @inbounds v.vec[ii]
    end

    @inline function Base.getindex(v::FractionVector{T,<:Real},j::Int) where T
        @boundscheck checkbounds(v, j)
        return ifelse(v.i==j,v.val,v.vec)
    end

    Base.IndexStyle(::Type{<:FractionVector}) = IndexLinear()

    function Solvers.primalval(v::FractionVector)
        vec₀ = Solvers.primalval(v.vec)
        val₀ = Solvers.primalval(v.val)
        return FractionVector(vec₀,val₀,v.i)
    end
end
