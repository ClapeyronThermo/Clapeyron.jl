module Compressed4DMatrices

using LinearAlgebra

"""
    Compressed4DMatrix{T}

A compressed storage that represents interactions between association sites of different components.

### Conceptual Model

Consider a system with `nc` components. Each component `i` has `n_i` association sites, numbered `1..n_i`. 
The full association matrix is of size `N × N`, where `N = sum(n_i)`, and its entries are indexed by `(component, site)` pairs:

    (i, a) ↔ (j, b)

The matrix is **symmetric** under exchange of the two pairs:

    A[(i,a), (j,b)] == A[(j,b), (i,a)]

To avoid storing duplicate entries, only one triangular half is kept:

- For `i < j`: store all `(a,b)` pairs (no symmetry within the block).
- For `i == j`: store only `a ≤ b` (self‑association symmetric).

The storage exploits the block structure: each component `i` has a contiguous range of site indices determined by `site_offsets`. 
The data is stored as a flat vector `values` and a parallel vector `indices` of encoded integers that uniquely identify each stored `(i,j,a,b)` tuple. 

### External Interface

- `m[i,j]` returns an `AssocView` that behaves like a dense matrix of size `n_i × n_j` (or `n_j × n_i` if the block is stored transposed). 
Indexing this view returns the appropriate value, handling symmetry and transposition automatically.
- `dropzeros!(m)` removes entries whose value is zero (the structural shape remains).

This type is used internally in Clapeyron.jl as the storage type for `AssocParam`.

## Constructors

Several constructors are provided to create a `Compressed4DMatrix` from different
input forms:

1. **From block sizes** – the most basic constructor:
   ```julia
   Compressed4DMatrix{T}(bsizes::AbstractVector{Int})
   ```
   where `bsizes[i]` is the number of association sites for component `i`.
   Creates an empty matrix (all values zero) with the specified block structure.

2. **From a matrix of matrices** – useful when a full block matrix is already assembled:
   ```julia
   Compressed4DMatrix(x::AbstractMatrix{<:AbstractMatrix{T}})
   ```
   `x` must be square (`nc × nc`), where `nc` is the number of components. 
   Each `x[i,j]` is a matrix of size `n_i × n_j` (or `n_j × n_i` if only one triangular half is provided; the constructor will handle transposition).
   The diagonal blocks `x[i,i]` must be square and symmetric (only the upper triangle is used).

3. **From explicit list of 4‑tuples** – to set specific entries:
   ```
   Compressed4DMatrix(vals::AbstractVector, ijab::AbstractVector{NTuple{4,<:Integer}})
   ```
   where each tuple `(i, j, a, b)` specifies a component‑site pair, and `vals` gives the corresponding value. 
   The entries can be in any order; they will be canonicalised and sorted internally. 
   Duplicate indices will overwrite the earlier value.

4. **From separate `ij` and `ab` vectors** – for convenience:
   ```
   Compressed4DMatrix(vals::AbstractVector, ij::Vector{NTuple{2,Int}}, ab::Vector{NTuple{2,Int}})
   ```
   similar to the constructor from `(i, j, a, b)` indices, but the component and site indices are in different vectors.

5. **Empty constructor**:
   ```
   Compressed4DMatrix{T}()
   ```
   creates an empty matrix with no components (useful as a placeholder).
"""
struct Compressed4DMatrix{T,V<:AbstractVector{T}}
    values::V #list of values 1:w, plain vector
    indices::Vector{Int} #global indices, 16x4
    site_offsets::Vector{Int} #start and end position of each block, n + 1, last index is nk, the size of the sum of all blocks.
end

Compressed4DMatrix(values::AbstractVector{T},indices::Vector{Int},site_offsets::Vector{Int}) where {T} = Compressed4DMatrix{T,typeof(values)}(values,indices,site_offsets)

#indexing:
#=
each i,j,a,b is stored in a packed format (16 bits each).

Having one indices vector helps simplifying the design.

There is also a site_offsets vector that enforces the shape, so bounds checking can now be performed.

The symmetry is enforced at construction time
- i != j, then i > j, no restrictions on a,b
- i == j, then a > b
=#

export Compressed4DMatrix

Base.length(m::Compressed4DMatrix) = length(m.indices)
Base.getindex(m::Compressed4DMatrix,idx::Int) = m.values[idx]
Base.eltype(m::Compressed4DMatrix{T}) where T = T
Base.setindex!(m::Compressed4DMatrix,val,i::Int) = Base.setindex!(m.values,val,i)

function Base.copyto!(dest::Compressed4DMatrix,src::Base.Broadcast.Broadcasted) #general, just copies the values, used in a .= f.(a)
    Base.copyto!(dest.values,src)
    return dest
end

function Base.copyto!(dest::Compressed4DMatrix,src::AbstractArray) #general, just copies the values, used in a .= f.(a)
    Base.copyto!(dest.values,src)
    return dest
end

function Base.copyto!(dest::Compressed4DMatrix,src::Compressed4DMatrix) #specific
    n = length(src.values)
    copyto!(resize!(dest.values,n),src.values)
    copyto!(resize!(dest.indices,n),src.indices)
    copyto!(resize!(dest.site_offsets,n),src.site_offsets)
    return dest
end

#also the number of components.
@inline nblocks(m::Compressed4DMatrix) = length(m.site_offsets) - 1

#size of the association matrix generated by the list of association site-pairs.
@inline matsize(m::Compressed4DMatrix) = last(m.site_offsets)

@inline blocksize(m::Compressed4DMatrix,i) = blocksize(m.site_offsets,i)

@inline function blocksize(offs::AbstractVector{T},i) where T <: Integer
    @boundscheck checkbounds(offs,i+1)
    @inbounds begin
        n0 = offs[i]
        n1 = offs[i+1]
        dn = n1 - n0
    end
end

#idx (1:length(m.values)) -> i,j,a,b
@inline function idx_to_ijab(m::Compressed4DMatrix,idx)
    ids = m.indices
    @boundscheck checkbounds(ids,idx)
    @inbounds begin
        return idx_to_ijab(ids[idx])
    end
end

@inline function idx_to_ijab(ijab::Union{Int64,UInt64})
    bjai_uint = reinterpret(NTuple{4,UInt16},ijab)
    b,a,j,i = map(Int,bjai_uint)
    return (i,j,a,b)
end


@inline function idx_to_ijab(ijab::Union{Int32,UInt32})
    bjai_uint = reinterpret(NTuple{4,UInt8},ijab)
    b,a,j,i = map(Int32,bjai_uint)
    return (i,j,a,b)
end

function ijab_to_idx(i::T,j::T,a::T,b::T) where T<:Union{Int64,UInt64}
    bjai_uint = map(UInt16,(b,a,j,i))
    return reinterpret(Int64,bjai_uint)
end

ijab_to_idx(ijab::NTuple{4,T}) where T = ijab_to_idx(ijab...)

@inline function canonical_index(ijab)
    i,j,a,b = ijab
    return canonical_index(i,j,a,b)
end

@inline function canonical_index(ij,ab)
    i,j = ij
    a,b = ab
    return canonical_index(i,j,a,b)
end

@inline function canonical_index(_i,_j,_a,_b)
    is_symmetric = _i == _j
    is_transpose = _i > _j
    i,j   = minmax(_i,_j)
    a,b = ifelse(is_symmetric,minmax(_a,_b),ifelse(_i != i,(_b,_a),(_a,_b)))
    return i,j,a,b
end

function IJABIterator(ik)
    idx,w = ik
    i,j,a,b = idx_to_ijab(w)
    return idx,(i,j),(a,b)
end

function indices(m::Compressed4DMatrix)
    return Iterators.map(IJABIterator,enumerate(m.indices))
end

function offsets_from_bsizes!(bsizes)
    nc = length(bsizes)
    resize!(bsizes, nc + 1)

    @inbounds for i in nc:-1:1
        bsizes[i+1] = bsizes[i]
    end
    
    bsizes[1] = 1
    offset = 1
    @inbounds for i in 1:nc
        bsize_i =  bsizes[i + 1]
        offset += bsize_i
        bsizes[i+1] = offset
    end
    return bsizes
end

function bsizes_from_offsets!(site_offsets::Vector{Int})
    nc = length(site_offsets) - 1

    # recover bsizes[i] = site_offsets[i+1] - site_offsets[i] + 1
    @inbounds for i in 1:nc
        site_offsets[i] = site_offsets[i+1] - site_offsets[i] + 1
    end

    resize!(site_offsets, nc)
    return site_offsets
end


#Constructs a ``Compressed4DMatrix` from a list of nc site sizes.

function Compressed4DMatrix{T}(bsizes::AbstractVector{Int}) where T
    nc = length(bsizes)
    #this field is exacly the same as the one stored in SiteParam.n_sites.p
    site_offsets = offsets_from_bsizes!(copy(bsizes))
    return c4d_from_site_offsets(T,site_offsets)
end

Compressed4DMatrix(bsizes::AbstractVector{Int}) = Compressed4DMatrix{Float64}(bsizes)

c4d_from_site_offsets(site_offsets) = c4d_from_site_offsets(Float64,site_offsets)
function c4d_from_site_offsets(::Type{T},site_offsets) where T
    nc = length(site_offsets) - 1
    nk = site_offsets[end]
    C = nk + 1  #compression base
    indices = Int[]

    for i in 1:nc
        ni = site_offsets[i]
        ni1 = site_offsets[i+1]
        dni = ni1 - ni#blocksize(m,i)
        iszero(dni) && continue
        for j in i:nc
            nj = site_offsets[j]
            nj1 = site_offsets[j+1]
            dnj = nj1 - nj #blocksize(m,j)
            iszero(dnj) && continue
            for a in 1:dni
                if i == j
                    _na = a
                else
                    _na = 1
                end
                for b in _na:dnj
                    push!(indices,ijab_to_idx(i,j,a,b))
                end
            end
        end
    end

    nv = length(indices)
    values = zeros(Float64,nv)
    return Compressed4DMatrix(values, indices, site_offsets)
end

function Base.show(io::IO,mime::MIME"text/plain",m::Compressed4DMatrix{T}) where T
    nv = length(m.values)
    print(io,typeof(m)," with ",nv," entr",(nv == 1 ? "y:" : "ies"))
    !iszero(nv) && println(io,":")
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

function Base.:(==)(p1::Compressed4DMatrix,p2::Compressed4DMatrix)
    return (p1.values == p2.values) & (p1.indices == p2.indices) && (p1.site_offsets == p2.site_offsets)
end

function Compressed4DMatrix{T}() where T
    return Compressed4DMatrix(T[],Int[],[0])
end

function Compressed4DMatrix(x::AbstractMatrix{<:AbstractMatrix{T}}) where T

    if length(x) == 0
        return Compressed4DMatrix{T}()
    end

    nc = LinearAlgebra.checksquare(x)
    iszero(nc) && return Compressed4DMatrix{T}()

    bsizes = zeros(Int,nc)
    #get the size of the association matrices from the pure componeent (diagonal) data.
    for i in 1:nc
        xii = x[i,i]
        ni = LinearAlgebra.checksquare(xii)
        bsizes[i] = ni
    end

    site_offsets = offsets_from_bsizes!(bsizes)
    mat = c4d_from_site_offsets(T,site_offsets)
    for (idx,(i,j),(a,b)) in indices(mat)
        xij = x[i,j]
        if iszero(prod(size(xij))) && i != j #off-diagonal was stored asymetrically, so we only use the non-empty side, whatever that is.
            xij2 = transpose(x[j,i])
            mat.values[idx] = xij2[b,a]
        elseif iszero(prod(size(xij))) && i == j
            #that means that we stored an empty index, error out
            error("invalid index (i,j,a,b) = $((i,j,a,b)) inside a CompressedAssocMatrix.")
        else
            mat.values[idx] = xij[a,b]
        end
    end
    return mat
end

function dropzeros!(mat::Compressed4DMatrix)
    #note, while indices were dropped, the actual structure, encoded in site_offsets, is kept intact
    nonzero_idx = findall(!iszero,mat.values)
    keepat!(mat.values,nonzero_idx)
    keepat!(mat.indices,nonzero_idx)
    return mat
end

struct AssocView{T,V<:Compressed4DMatrix{T}} <: AbstractMatrix{T}
    m::V
    ijab::Int
    C::Int
    function AssocView(::Nothing,values::Compressed4DMatrix{T},ijab::Int,C::Int) where {T}
        return new{T,typeof(values)}(values,ijab,C)
    end
end

@inline function AssocView(m::Compressed4DMatrix,i,j)
    nk = matsize(m)
    @boundscheck begin
        checkbounds(Base.OneTo(nk),i)
        checkbounds(Base.OneTo(nk),j)
    end
    @inbounds begin
        C = nk
        a,b = blocksize(m,i),blocksize(m,j)
        return AssocView(nothing,m,ijab_to_idx(i,j,a,b),C)
    end
end

@inline function Base.getindex(m::Compressed4DMatrix,i,j)
    nk = matsize(m)
    @boundscheck begin
        checkbounds(Base.OneTo(nk),i)
        checkbounds(Base.OneTo(nk),j)
    end

    @inbounds begin
        return AssocView(m,i,j)
    end
end



@inline function Base.size(m::AssocView)
    i,j,a,b = idx_to_ijab(m.ijab)
    return ifelse(i >= j,(a,b),(b,a))
end

function Base.summary(io::IO,m::AssocView{T}) where T
    s1,s2 = size(m)
    print(io,s1,"×",s2)
    print(io," ")
    i,j,a,b = idx_to_ijab(m.ijab)
    if i == j
        print(io,"symmetric ")
    end

    if i > j
        print(io,"transposed")
    end
    print(io,typeof(m))
end

function Base.print_array(io::IO, m::AssocView{T}) where T
    s = size(m)
    M = Matrix{T}(undef,s)
    M .= 0
    for i in 1:s[1]
        for j in 1:s[2]
            M[i] = m[i,j]
        end
    end
    Base.print_array(io,M)
end


Base.eltype(m::AssocView{T}) where T = T

#returns the absolute index. that is. it is directly indexable by the parent array
@inline function validindex(m::AssocView{TT},k1::Int,k2::Int) where TT
    ijab = m.ijab
    _i,_j,_a,_b = idx_to_ijab(ijab) #unwrap the ijab, where ij is the specified ij, and ab is the size of the block ij
    (iszero(_a) || iszero(_b)) && return zero(eltype(w)) #zero--sized assoc view has no indices
    if _i > _j
        k1,k2 = k2,k1 #swap indices if we are in the presence of a transpose position
    end
    idxs = m.m.indices
    i,j,si,sj = canonical_index(_i,_j,k1,k2) #canonical index: returns the form that we store.
    k = ijab_to_idx(i,j,si,sj) #transform to compressed index
    len = length(idxs)
    T = eltype(idxs)
    w = searchsortedfirst(idxs,k) #search unique compressed index
    w > len && return T(0)
    in_idxs = @inbounds(idxs[w]) == k
    return ifelse(in_idxs,T(w),T(0))
end

@inline function validindex(m::Compressed4DMatrix{TT},ijab::NTuple{4,Int}) where TT
    idxs = m.indices
    can = canonical_index(ijab)
    k  = ijab_to_idx(canonical_index(ijab))
    len = length(idxs)
    T = eltype(idxs)
    w = searchsortedfirst(idxs,k,T(1),T(len),Base.Order.ForwardOrdering()) #search unique compressed index
    w > len && return T(0)
    in_idxs = @inbounds(idxs[w]) == k
    return ifelse(in_idxs,T(w),T(0))
end

@inline function Base.getindex(m::AssocView{T},i::Int,j::Int) where T
    idx = validindex(m,i,j)
    iszero(idx) && return zero(T)
    @inbounds begin
        return m.m.values[idx]
    end
end

@inline function Base.setindex!(m::AssocView{T},value,a::Int,b::Int) where T
    idx = validindex(m,a,b)
    vals = m.m.values
    @boundscheck begin
        checkbounds(m.m.indices,idx)
    end
    @inbounds begin 
        vals[idx] = value
    end
end

@inline function Base.setindex!(m::Compressed4DMatrix{T},value,ijab::NTuple{4,Int}) where T
    idx = validindex(m,ijab)
    vals = m.values
    @boundscheck begin
        checkbounds(m.indices,idx)
    end
    @inbounds begin 
        vals[idx] = value
    end
end


@inline function Base.getindex(m::Compressed4DMatrix{T},ijab::NTuple{4,Int}) where T
    idx = validindex(m,ijab)
    vals = m.values
    iszero(idx) && return zero(T)
    @boundscheck begin
        checkbounds(m.indices,idx)
    end
    @inbounds begin 
        return vals[idx]
    end
end

function site_indices(m::Compressed4DMatrix, w::Int)
    i,j,a,b = idx_to_ijab(m, w)   # decode (i,j,a,b)
    # global site index for component i, site a
    ia = m.site_offsets[i] + a - 1
    jb = m.site_offsets[j] + b - 1
    return ia, jb
end

function Base.Matrix(m::Compressed4DMatrix{T}) where T
    N = matsize(m)   # total number of sites
    A = zeros(T, N, N)
    vals = m.values
    for (w, val) in enumerate(vals)
        @inbounds begin
            ia, jb = _site_indices(m, w)
            A[ia, jb] = val
            A[jb, ia] = val # symmetry
        end   
    end
    return A
end

#compatibility functions

function infer_site_offsets(ij::AbstractVector,ab::AbstractVector,nc = 0)
    if iszero(nc)
        _nc = maximum(max,ij)
    else
        _nc = nc
        bsizes = zeros(Int,nc)
        for k in 1:length(ij)
            i,j = ij[k]
            a,b = ab[k]
            bsizes[i] = max(bsizes[i],a)
            bsizes[j] = max(bsizes[j],b)
        end
    end
    return offsets_from_bsizes!(bsizes)
end

function infer_site_offsets(ijab,nc = 0)
    if iszero(nc)
        _nc = maximum(max,ij)
    else
        _nc = nc
        bsizes = zeros(Int,nc)
        for k in 1:length(ijab)
            i,j,a,b = ijab[k]
            bsizes[i] = max(bsizes[i],a)
            bsizes[j] = max(bsizes[j],b)
        end
    end
    return offsets_from_bsizes!(bsizes)
end


function Compressed4DMatrix(vals::AbstractVector{T},ijab::AbstractVector{NTuple{4,Int}},offs::AbstractVector{Int} = infer_site_offsets(ijab)) where {T}
    m = c4d_from_site_offsets(offs)
    for (idx,_ijab) in pairs(ijab)
        i,j,a,b = canonical_index(_ijab)
        m[(i,j,a,b)] = vals[idx]
    end
    dropzeros!(m)
    return m
end

function Compressed4DMatrix(vals::AbstractVector{T},ij::Vector{NTuple{2,Int}},ab::Vector{NTuple{2,Int}},unsafe::Bool) where T
    if unsafe
        @warn "Compressed4DMatrix: unsafe keyword is deprecated."
    end
    return Compressed4DMatrix(vals,ij,ab)
end

function Compressed4DMatrix(vals::AbstractVector{T},ij::Vector{NTuple{2,Int}},ab::Vector{NTuple{2,Int}},offs::AbstractVector{Int} = infer_site_offsets(ijab)) where T
    m = Compressed4DMatrix{T}(offs)
    for (idx,(_ij,_ab)) in pairs(zip(ij,ab))
        i,j,a,b = canonical_index(_ij,_ab)
        @inbounds begin
            m[i,j][a,b] = vals[idx]
        end
    end
    dropzeros!(m)
    return m
end

assoc_is_pure(m::Compressed4DMatrix,i) = assoc_is_pure(m.indices[i])

function assoc_is_pure(ix)
    i,j,a,b = idx_to_ijab(ix)
    return i == j
end

assoc_self_assoc(m::Compressed4DMatrix,i) = assoc_is_pure(m.indices[i])

function assoc_self_assoc(ix)
    i,j,a,b = idx_to_ijab(ix)
    return i == j && a == b
end

end #module

using .Compressed4DMatrices
using .Compressed4DMatrices: AssocView, indices, dropzeros!, validindex, idx_to_ijab, ijab_to_idx, assoc_is_pure, assoc_self_assoc
import .Compressed4DMatrices: Compressed4DMatrix

"""
    assoc_similar(mat::Compressed4DMatrix)
    assoc_similar(mat::Compressed4DMatrix,::Type{𝕋}) where 𝕋 <:Number
    
Returns a `Clapeyron.Compressed4DMatrix` of the same shape as the input, with the same element type as `𝕋`.
"""
function assoc_similar(m::Compressed4DMatrix,::Type{𝕋}) where 𝕋 <:Number
    newvalues = zeros(𝕋,length(m.values))
    return Compressed4DMatrix(newvalues,m.indices,m.site_offsets)
end

assoc_similar(mat::Compressed4DMatrix{T}) where T = assoc_similar(mat,T)

function Solvers.primalval(x::Compressed4DMatrix{T}) where T
    return Compressed4DMatrix(Solvers.primalval(x.values),x.indices,x.site_offsets)
end

function Solvers.primalval_eager(x::Compressed4DMatrix{T}) where T
    return Compressed4DMatrix(Solvers.primalval_eager(x.values),x.indices,x.site_offsets)
end

#Y = A*X + Y
function LinearAlgebra.mul!(Y::AbstractVector, A::Compressed4DMatrix, X::AbstractVector, α::Number, β::Number)
    N = matsize(A)
    @assert length(Y) == N && length(X) == N
    Y .*= β
    @inbounds for (idx, (i, j), (a, b)) in indices(A)
        v = A.values[idx] * α
        ia = A.site_offsets[i] + a - 1
        jb = A.site_offsets[j] + b - 1
        Y[ia] += v * X[jb]
        if !(i == j && a == b)   # avoid double‑counting diagonal self‑term
            Y[jb] += v * X[ia]
        end
    end
    return Y
end
