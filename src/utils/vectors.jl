struct FractionVector{T,V} <: AbstractVector{T}
    vec::V
    val::T
end
#=

Fraction Vector
useful when expressing fractions in n-1 components.
the last component is calculated at build time.
it allocates less than creating a new vector or appending.
=#
##
function FractionVector(v::AbstractVector)
    a = one(eltype(v))
    # any(x->x<0,v) && throw(DomainError(v,"all elements of a fraction vector should be positive."))
    a -=sum(v)
    # a < 0 && throw(DomainError(a,"the values of the input vector add to more than one"))
    return FractionVector(v,a) 
end

function FractionVector(v::Real)
    a = one(v) 
    # (v < zero(v)) && throw(DomainError(v,"all elements of a fraction vector should be positive."))
    a -= v
    # a < 0 && throw(DomainError(a,"the values of the input vector add to more than one"))
    return FractionVector(v,a)
end


@inline Base.eltype(v::FractionVector{T}) where T = T
@inline Base.length(v::FractionVector)::Int = Base.length(v.vec) + 1

@inline function Base.length(v::FractionVector{T,<:Real})::Int where T
    return 2
end

@inline Base.size(v::FractionVector) = (length(v),)
@inline function Base.getindex(v::FractionVector,i::Int)
    @boundscheck checkbounds(v, i)
    return ifelse(length(v)==i,v.val,v.vec[min(length(v.vec),i)])
end

@inline function Base.getindex(v::FractionVector{T,<:Real},i::Int) where T
    @boundscheck checkbounds(v, i)
    return ifelse(i==1,v.vec,v.val)
end

Base.IndexStyle(::Type{<:FractionVector}) = IndexLinear()

##

struct CompressedAssocMatrix{T} 
    values::Vector{T}
    outer_indices::Vector{Tuple{Int,Int}} #index of components
    inner_indices::Vector{Tuple{Int,Int}} #index of sites
    outer_size::Tuple{Int,Int} #size of component matrix
    inner_size::Tuple{Int,Int} #size of sites matrices
end

function Base.show(io::IO,mime::MIME"text/plain",m::CompressedAssocMatrix{T}) where T
    n = length(m.values)
    println(io,typeof(m)," with ",n," entr",(n == 1 ? "y:" : "ies:"))
    for (idx,(i,j),(a,b)) in indices(m)
        print(io," ",(i,a)," >=< ",(j,b),": ",m.values[idx])
    end
end

function Base.show(io::IO,m::CompressedAssocMatrix{T}) where T
    print(io,typeof(m))
    print(io,m.values)
end
function CompressedAssocMatrix{T}() where T
    return CompressedAssocMatrix(T[],Tuple{Int,Int}[],Tuple{Int,Int}[],(0,0),(0,0))
end

const MatrixofMatrices{T} = AbstractMatrix{<:AbstractMatrix{T}} where T

function CompressedAssocMatrix(x::MatrixofMatrices{T}) where T
    outer_size = size(x)
    is1, is2 = 0, 0
    os1, os2 = outer_size
    for i in eachindex(x)
        _is1,_is2 = size(x[i])
        is1,is2 = max(_is1,is1),max(_is2,is2)
    end
    inner_size = (is1,is2)

    values = T[]
    inner_indices = Tuple{Int,Int}[]
    outer_indices = Tuple{Int,Int}[]

    if iszero(os1) & iszero(os2)
        return CompressedAssocMatrix(values,outer_indices,inner_indices,outer_size,inner_size)
    end
    

    for i in 1:os1
        for j in i:os2 #there can be self association
            xi = x[i,j]
            if iszero(length(xi))
                continue
            end
            a1,a2 = size(xi)
            for a in 1:a1
                start = ifelse(i == j,a,1)
                for b in start:a2 #but not on the same site
                    if !_iszero(xi[a,b]) 
                        push!(values,xi[a,b])
                        push!(outer_indices,(i,j))
                        push!(inner_indices,(a,b))
                    end
                end
            end
        end
    end
    return CompressedAssocMatrix(values,outer_indices,inner_indices,outer_size,inner_size)
end

function Base.getindex(m::CompressedAssocMatrix,i::Int,j::Int)
    i,j = minmax(i,j)
    @inbounds begin
    idx = searchsorted(m.outer_indices,(i,j))
    return AssocView(view(m.values,idx),view(m.inner_indices,idx),m.inner_size)
    end
end

#Base.eltype(m::CompressedAssocMatrix{T}) where T = T

function Base.setindex!(m::CompressedAssocMatrix,val,i::Int)
    @inbounds begin
        m.values[i] = val
    end
end

struct AssocView{T,V,I} <: AbstractMatrix{T}
    values::V
    indices::I
    size::Tuple{Int,Int}
end

function AssocView(V::Vi,I::Ii,s) where Ii where Vi <:AbstractVector{<:T} where T
    return AssocView{T,Vi,Ii}(V,I,s)
end

Base.size(m::AssocView) = m.size
Base.eltype(m::AssocView{T}) where T = T

function Base.getindex(m::AssocView{T},i::Int,j::Int) where T
    @inbounds begin
        idxs = searchsorted(m.indices,(i,j))
        if iszero(length(idxs)) 
            idxs = searchsorted(m.indices,(j,i))
            iszero(length(idxs)) &&  return zero(T)  
        end
        idx = first(idxs)
        return m.values[idx]
    end
end

function zero_assoc(m::CompressedAssocMatrix,::Type{T} = Float64) where T <:Number
    newvalues = zeros(T,length(m.values))
    return CompressedAssocMatrix(newvalues,m.outer_indices,m.inner_indices,m.outer_size,m.inner_size)
end

function zero_assoc(m::Matrix{<:Matrix},a::Type{T} = Float64) where T <:Number
    return zero_assoc(CompressedAssocMatrix(m),a)
end

#=
function zero_assoc(m::CompressedAssocMatrix,x::T = 0.0) where T <:Number
    newvalues = fill(x,length(m.values))
    return CompressedAssocMatrix(newvalues,m.outer_indices,m.inner_indices,m.outer_size,m.inner_size)
end
=#
function indices(x::CompressedAssocMatrix)
    return zip(1:length(x.values),x.outer_indices,x.inner_indices)
end

function indices(x::PackedVofV)
    return x.p
end