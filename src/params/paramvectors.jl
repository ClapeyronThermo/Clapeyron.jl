#=
The Vectors defined here are necessary to 
ClapeyronParam
=#
struct Compressed4DMatrix{T,V<:AbstractVector{T}} 
    values::V
    outer_indices::Vector{Tuple{Int,Int}} #index of components
    inner_indices::Vector{Tuple{Int,Int}} #index of sites
    outer_size::Tuple{Int,Int} #size of component matrix
    inner_size::Tuple{Int,Int} #size of sites matrices
end

function Base.show(io::IO,mime::MIME"text/plain",m::Compressed4DMatrix{T}) where T
    n = length(m.values)
    println(io,typeof(m)," with ",n," entr",(n == 1 ? "y:" : "ies:"))
    for (idx,(i,j),(a,b)) in indices(m)
        if idx != 1
        println(io)
        end
        print(io," ",(i,a)," >=< ",(j,b),": ",m.values[idx])
    end
end

function Base.show(io::IO,m::Compressed4DMatrix{T}) where T
    print(io,typeof(m))
    print(io,m.values)
end
function Compressed4DMatrix{T}() where T
    return Compressed4DMatrix(T[],Tuple{Int,Int}[],Tuple{Int,Int}[],(0,0),(0,0))
end

const MatrixofMatrices{T} = AbstractMatrix{<:AbstractMatrix{T}} where T

function Compressed4DMatrix(x::MatrixofMatrices{T}) where T
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
        return Compressed4DMatrix(values,outer_indices,inner_indices,outer_size,inner_size)
    end

    #self association

    for i in 1:os1
        for j in i:os2 #there can be self association
            xi = x[i,j]
            if iszero(length(xi))
                continue
            end
            a1,a2 = size(xi)
            for a in 1:a1
                start = ifelse(i == j,a,1)
                for b in start:a2 #this includes (i,i)(a,a) (sCKSAFT)
                    if !_iszero(xi[a,b]) 
                        push!(values,xi[a,b])
                        push!(outer_indices,(i,j))
                        push!(inner_indices,(a,b))
                    end
                end
            end
        end
    end
    return Compressed4DMatrix{T,Vector{T}}(values,outer_indices,inner_indices,outer_size,inner_size)
end

function Base.getindex(m::Compressed4DMatrix,i::Int,j::Int)
    i,j = minmax(i,j)
    @inbounds begin
    idx = searchsorted(m.outer_indices,(i,j))
    return AssocView(view(m.values,idx),view(m.inner_indices,idx),m.inner_size)
    end
end

#Base.eltype(m::Compressed4DMatrix{T}) where T = T

function Base.setindex!(m::Compressed4DMatrix,val,i::Int)
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

function zero_assoc(m::Compressed4DMatrix,::Type{T} = Float64) where T <:Number
    newvalues = zeros(T,length(m.values))
    return Compressed4DMatrix(newvalues,m.outer_indices,m.inner_indices,m.outer_size,m.inner_size)
end

function zero_assoc(m::Matrix{<:Matrix},a::Type{T} = Float64) where T <:Number
    return zero_assoc(Compressed4DMatrix(m),a)
end

#=
function zero_assoc(m::Compressed4DMatrix,x::T = 0.0) where T <:Number
    newvalues = fill(x,length(m.values))
    return Compressed4DMatrix(newvalues,m.outer_indices,m.inner_indices,m.outer_size,m.inner_size)
end
=#
function indices(x::Compressed4DMatrix)
    xin = x.outer_indices
    l = 1:length(xin) 
    return zip(l,xin,x.inner_indices)
end

function indices(x::PackedVofV)
    return x.p
end


struct SparsePackedMofV{E,P<:PackedVofV}<:SparseArrays.AbstractSparseMatrixCSC{E,Int}
    storage::P
    idx::SparseMatrixCSC{Int,Int}
end

function SparsePackedMofV(storage,idx)
    E = eltype(storage)
    P = typeof(storage)
    return SparsePackedMofV{E,P}(storage,idx)
end

_findnz(x::SparseMatrixCSC) = SparseArrays.findnz(x)
function _findnz(x::AbstractMatrix{<:AbstractVector})
    idxs = findall(z->!iszero(length(z)),x)
    i = first.(idxs)
    j = last.(idxs)
    vals = x[idxs]
    return i,j,vals
end

function SparsePackedMofV(m_of_v::AbstractMatrix{<:AbstractVector})
    i,j,unpack_sparse_vals = _findnz(m_of_v)
    #unpack_sparse_vals = m_of_v[sparse_idx]
    pack_sparse_vals = PackedVectorsOfVectors.pack(unpack_sparse_vals)
    len_linear = length(unpack_sparse_vals)
    idx_linear = collect(1:len_linear)
    m,n = size(m_of_v)
    idx = sparse(i,j,idx_linear,m,n)
    E = eltype(pack_sparse_vals)
    P = typeof(pack_sparse_vals)
    return SparsePackedMofV{E,P}(pack_sparse_vals,idx)
end

function Base.getindex(x::SparsePackedMofV,i::Int,j::Int)
    _idx = x.idx[i,j]
    iszero(_idx) && (return view(x.storage.v,1:0))
    
    return x.storage[_idx]
end
@inline Base.size(x::SparsePackedMofV) = size(x.idx)
@inline SparseArrays.nnz(x::SparsePackedMofV) = length(x.storage.p)
@inline function SparseArrays.findnz(x::SparsePackedMofV)
    i,j,_ = SparseArrays.findnz(x.idx)
    return i,j,x.storage
end
@inline SparseArrays.nonzeros(x::SparsePackedMofV) = x.storage
@inline SparseArrays.rowvals(x::SparsePackedMofV) = SparseArrays.rowvals(x.idx)
@inline SparseArrays.nzrange(A::SparsePackedMofV, col::Integer) = SparseArrays.nzrange(A.idx,col)
@inline SparseArrays.getcolptr(A::SparsePackedMofV) = SparseArrays.getcolptr(A.idx)
export SparsePackedMofV

function Base.show(io::IO,::MIME"text/plain",A::SparsePackedMofV)
    m,n = size(A)
    println(io,"$(m)×$(n) Sparse Packed Matrix of Vectors of eltype $(eltype(A.storage)) with $(length(A.storage)) non empty values:")
    vals = A.storage
    rows = rowvals(A)
    for j in 1:n
        for ii ∈ nzrange(A, j)
            i = rows[ii]
            val= vals[ii]
            println(io,"  ($i,$j) => $val")
        end
    end
end
