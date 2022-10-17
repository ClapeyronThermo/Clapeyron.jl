"""
    Compressed4DMatrix{T,V<:AbstractVector{T}}

Struct used to hold association data. as its name says, it is a compressed 4D matrix containing all the non-zero combinations of component-site pairs.
The component-site pairs `(i,j,a,b)` are sorted lexicographically. the `(i,j)` pairs are stored in the `outer_indices` field, whereas the `(a,b)` pairs are stored in the `inner_indices` field. 

Let's see an associating model:
```julia-repl
julia> model = PCSAFT(["water","methanol","ethane"],assoc_options = AssocOptions(combining = :esd))
PCSAFT{BasicIdeal} with 3 components:
 "water"
 "methanol"
 "ethane"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol
```

We check out the `bondvol` parameter. note how ethane does not appear in the list:

```julia-repl
julia> model.params.bondvol
AssocParam{Float64}["water", "methanol", "ethane"]) with 4 values:
("water", "e") >=< ("water", "H"): 0.034868
("methanol", "e") >=< ("water", "H"): 0.03495053755004983
("methanol", "H") >=< ("water", "e"): 0.03495053755004983
("methanol", "e") >=< ("methanol", "H"): 0.035176
```

The underlying structure used to store `AssocParam` values is a `Compressed4DMatrix`:
```julia-repl
julia> vals = model.params.bondvol.values
Clapeyron.Compressed4DMatrix{Float64, Vector{Float64}} with 4 entries:
 (1, 1) >=< (1, 2): 0.034868
 (2, 1) >=< (1, 2): 0.03495053755004983
 (2, 2) >=< (1, 1): 0.03495053755004983
 (2, 1) >=< (2, 2): 0.035176

julia> vals.values
4-element Vector{Float64}:
 0.034868
 0.03495053755004983
 0.03495053755004983
 0.035176

julia> vals.outer_indices
4-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (2, 1)
 (2, 1)
 (2, 2)

julia> vals.inner_indices
4-element Vector{Tuple{Int64, Int64}}:
 (1, 2)
 (1, 2)
 (2, 1)
 (1, 2)
```
If we check the indices:
```julia-repl
julia> idxs = [(ij...,ab...) for (ij,ab) in zip(vals.outer_indices,vals.inner_indices)]
4-element Vector{NTuple{4, Int64}}:
 (1, 1, 1, 2)
 (2, 1, 1, 2)
 (2, 1, 2, 1)
 (2, 2, 1, 2)

julia> issorted(idxs)
true
```
"""
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
    indices = Tuple{Int,Int,Int,Int}[]

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
                        push!(indices,(i,j,a,b))
                    end
                end
            end
        end
    end
    idx =  sortperm(indices)
    indices = indices[idx]
    outer_indices = [(c[1],c[2]) for c ∈ indices]
    inner_indices = [(c[3],c[4]) for c ∈ indices]
    values = values[idx]
    return Compressed4DMatrix{T,Vector{T}}(values,outer_indices,inner_indices,outer_size,inner_size)
end

"""
    Compressed4DMatrix(vals::AbstractVector,idxs::AbstractVector)

builds a `Compressed4DMatrix`, using a list of values and a list of `(i,j,a,b)::NTuple{4,Int}` indices:

## Example:
```
julia> ijab, vals = [(1,1,1,2)], [3.0]
([(1, 1, 1, 2)], [3.0])

julia> Clapeyron.Compressed4DMatrix(vals,ijab)
Clapeyron.Compressed4DMatrix{Float64, Vector{Float64}} with 1 entry:
 (1, 1) >=< (1, 2): 3.0
```
"""
function Compressed4DMatrix(vals::AbstractVector,idxs::AbstractVector)
    if issorted(idxs)
        new_vals = copy(vals)
        new_idxs = copy(idxs)
    else
        sort_idxs = sortperm(idxs)
        new_idxs = idxs[sort_idxs]
        new_vals = vals[sort_idxs]
    end
    ij = map(x -> (x[1],x[2]),new_idxs)
    ab = map(x -> (x[3],x[4]),new_idxs)
    return Compressed4DMatrix(new_vals,ij,ab,true)
end

"""
    Compressed4DMatrix(vals,ij,ab,unsafe::Bool = false)

Builds a `Compressed4DMatrix`, using a list of values, a list of `ij:Tuple{Int,Int}` outer indices and a list of `ab:Tuple{Int,Int}` inner indices.

If `unsafe` is true, `ij` and `ab` will be considered sorted, and will build a `Compressed4DMatrix` directly, using the same reference to `vals`, `ij` and `ab`:
## Example:
```julia-repl
julia> ij, ab, vals = [(1,1)], [(1,2)], [3.0]
([(1, 1)], [(1, 2)], [3.0])

julia> assoc1,assoc2 = Clapeyron.Compressed4DMatrix(vals,ij,ab),Clapeyron.Compressed4DMatrix(vals,ij,ab,true)
(Clapeyron.Compressed4DMatrix{Float64, Vector{Float64}}[3.0], Clapeyron.Compressed4DMatrix{Float64, Vector{Float64}}[3.0])

julia> assoc1.values[1] = 100; (vals,assoc1.values[1])
([3.0], 100.0)

julia> assoc2.values[1] = 100; (vals,assoc2.values[1])
([100.0], 100.0)
```
"""
function Compressed4DMatrix(vals,ij,ab,unsafe::Bool = false)
    if !unsafe
        ijab = [(ij...,ab...) for (ij,ab) in zip(ij,ab)] 
        return Compressed4DMatrix(vals,ijab)
    end
    _ij_size = maximum(maximum(i) for i ∈ ij)
    _ab_size = maximum(maximum(i) for i ∈ ab)
    ij_size = (_ij_size,_ij_size)
    ab_size = (_ab_size,_ab_size)
    return Compressed4DMatrix(vals,ij,ab,ij_size,ab_size)
end
function Base.getindex(m::Compressed4DMatrix,i::Int,j::Int)
    j,i = minmax(i,j)
    @inbounds begin
    idx = searchsorted(m.outer_indices,(i,j))
    if iszero(idx)
        idx = searchsorted(m.outer_indices,(j,i))
    end
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

"""
    assoc_similar(mat::Compressed4DMatrix)
    assoc_similar(mat::Compressed4DMatrix,::Type{𝕋}) where 𝕋 <:Number)

returns a `Clapeyron.Compressed4DMatrix` of the same shape as the input, with the same element type as `𝕋`
"""
function assoc_similar(m::Compressed4DMatrix,::Type{𝕋}) where 𝕋 <:Number
    newvalues = similar(m.values,𝕋)
    return Compressed4DMatrix(newvalues,m.outer_indices,m.inner_indices,m.outer_size,m.inner_size)
end

assoc_similar(mat::Compressed4DMatrix{T}) where T = assoc_similar(mat,T)

function indices(x::Compressed4DMatrix)
    xin = x.outer_indices
    l = 1:length(xin) 
    return zip(l,xin,x.inner_indices)
end

"""
    SparsePackedMofV{T,V<:AbstractVector{T}} <:SparseArrays.AbstractSparseMatrixCSC{E,Int}

Sparse Matrix struct used internally to store a matrix of Vectors efficiently.

"""
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

#= Operations
function Base.:(+)(matrix::Compressed4DMatrix, x::Number)
    values = matrix.values .+ x
    return Compressed4DMatrix(values, matrix.outer_indices, matrix.inner_indices, matrix.outer_size, matrix.inner_size)
end

function Base.:(*)(matrix::Compressed4DMatrix, x::Number)
    values = matrix.values .* x
    return Compressed4DMatrix(values, matrix.outer_indices, matrix.inner_indices, matrix.outer_size, matrix.inner_size)
end

function Base.:(^)(matrix::Compressed4DMatrix, x::Number)
    values = matrix.values .^ x
    return Compressed4DMatrix(values, matrix.outer_indices, matrix.inner_indices, matrix.outer_size, matrix.inner_size)
end =#