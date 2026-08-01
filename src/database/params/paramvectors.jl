module Compressed4DMatrices

using Clapeyron: _iszero, _zero, parameterless_type
using LinearAlgebra

function low_color(text::AbstractString)
    colors = Base.text_colors
    g = colors[:light_black]
    reset = colors[:normal]
    return g * text * reset
end

low_color(symbol::Symbol) = low_color(":" * string(symbol))
low_color(x) = low_color(string(x))

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
    nv = length(src.values)
    ns = length(src.site_offsets)
    copyto!(resize!(dest.values,nv),src.values)
    copyto!(resize!(dest.indices,nv),src.indices)
    copyto!(resize!(dest.site_offsets,ns),src.site_offsets)
    return dest
end

#also the number of components.
@inline nblocks(m::Compressed4DMatrix) = length(m.site_offsets) - 1

#size of the association matrix generated by the list of association site-pairs.
@inline matsize(m::Compressed4DMatrix) = last(m.site_offsets) - 1

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

function insert_at!(x::Vector, index::Int, value)
    n = length(x)
    @boundscheck begin
        if index < 1 || index > n + 1
            throw(BoundsError(x, index))
        end
    end
    # Increase vector size by one
    resize!(x, n + 1)
    if !iszero(n) #handle n = 0
        xl = x[end-1]
        x[end] = xl
    end
    # Shift elements from the insertion point to the right
    for i in n:-1:(index+1)
        x[i] = x[i-1]
    end
    # Place the new value
    x[index] = value

    return x
end

function IJABIterator(ik)
    idx,w = ik
    i,j,a,b = idx_to_ijab(w)
    return idx,(i,j),(a,b)
end

indices(m::Compressed4DMatrix) = indices(m.indices)

function indices(m::AbstractVector{T}) where T <: Integer
    return Iterators.map(IJABIterator,enumerate(m))
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
        site_offsets[i] = site_offsets[i+1] - site_offsets[i]
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
    values = Vector{T}(undef,nv)
    _0 = _zero(T)
    values .= _0
    return Compressed4DMatrix(values, indices, site_offsets)
end

#if set_vals is set to false, then only the resizing and index creation is done, storing old values at the end
#if set_vals is set to true, then, additionally, the old values will be put in their corresponding positions of the extended matrix.
#in this way, we can "change" the index type completely via external functions, like what gc_to_comp_sites does. 
function __extend!(m::Compressed4DMatrix{T},set_vals::Bool) where T
    nv = length(m.indices)
    _set_vals = set_vals && nv != 0
    n = matsize(m)
    nf = div(n*(n+1),2)
    nc = nblocks(m)
    site_offsets = m.site_offsets
    ntt = nf + nv #extended size
    values = m.values
    indices = m.indices
    resize!(values,ntt)
    resize!(indices,ntt)

    #move values to the end.
    for i in 1:nv
        values[i + nf]  =  values[i]
        indices[i + nf] = indices[i]
    end
    w = nf + 1
    wf = 0
    _0 = _zero(T)

    #if we are setting vals, we expect a non-zero entries, with zero entries we just do not set values in the extended matrix.
    if _set_vals
        idx = indices[w]
        val = values[w]
    else
        idx = one(eltype(indices))
        val = _0
    end

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
                _na =  i == j ? a : 1
                for b in _na:dnj
                    wf += 1
                    ijab_idx = ijab_to_idx(i,j,a,b)
                    indices[wf] = ijab_idx
                    if ijab_idx == idx && w <= ntt && _set_vals
                        values[wf] = val
                        w += 1
                        idx = indices[min(ntt,w)]
                        val = values[min(ntt,w)]
                    else
                        values[wf] = _0
                    end
                end
            end
        end
    end
    return m,nf,nv
end

function extend!(m::Compressed4DMatrix)
    m,nf,nv = __extend!(m,true)
    resize!(m.values,nf)
    resize!(m.indices,nf)
    return m
end

function Base.show(io::IO,mime::MIME"text/plain",m::Compressed4DMatrix{T}) where T
    nv = length(m.values)
    print(io,typeof(m)," with ",nv," entr",(nv == 1 ? "y:" : "ies"))
    !iszero(nv) && println(io,":")
    is_str = T <: AbstractString
    for (idx,(i,j),(a,b)) in indices(m)
        if idx != 1
        println(io)
        end
        v = m.values[idx]
        print(io," ",(i,a)," >=< ",(j,b),": ")
        is_str ? print(io,"\"",v,"\"") : print(io,v)

    end
end

function Base.show(io::IO,m::Compressed4DMatrix{T}) where T
    print(io,parameterless_type(m))
    if get(io,:Clapeyron_eos_repr,false)
        print(io,'(')
        idxs = idx_to_ijab.(m.indices)
        print(io,m.values)
        print(io,',')
        print(io,idxs)
        print(io,',')
        print(io,m.site_offsets)
        print(io,')')
    else
        print(io,'{')
        print(io,T)
        print(io,'}')
        print(io,m.values)
    end
end

function Base.:(==)(p1::Compressed4DMatrix,p2::Compressed4DMatrix)
    return (p1.values == p2.values) & (p1.indices == p2.indices) && (p1.site_offsets == p2.site_offsets)
end

function Compressed4DMatrix{T}() where T
    return Compressed4DMatrix(T[],Int[],[0])
end

function Compressed4DMatrix(x::AbstractMatrix{<:AbstractMatrix{T}}) where T

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
            mat.values[idx] = xij2[a,b]
        elseif iszero(prod(size(xij))) && i == j
            #that means that we stored an empty index, error out
            throw(error("invalid index (i,j,a,b) = $((i,j,a,b)) inside a CompressedAssocMatrix."))
        else
            mat.values[idx] = xij[a,b]
        end
    end
    dropzeros!(mat)
    return mat
end

function dropzeros!(mat::Compressed4DMatrix)
    #note, while indices were dropped, the actual structure, encoded in site_offsets, is kept intact
    nonzero_idx = findall(!_iszero,mat.values)
    keepat!(mat.values,nonzero_idx)
    keepat!(mat.indices,nonzero_idx)
    return mat
end

struct AssocView{T,V<:Compressed4DMatrix{T}} <: AbstractMatrix{T}
    m::V
    ijab::Int
    function AssocView(values::Compressed4DMatrix{T},ijab::Int) where {T}
        return new{T,typeof(values)}(values,ijab)
    end
end

@inline function AssocView(m::Compressed4DMatrix,i,j)
    nc = nblocks(m)
    @boundscheck begin
        checkbounds(Base.OneTo(nc),i)
        checkbounds(Base.OneTo(nc),j)
    end
    @inbounds begin
        a,b = blocksize(m,i),blocksize(m,j)
        return AssocView(m,ijab_to_idx(i,j,a,b))
    end
end

@inline function Base.getindex(m::Compressed4DMatrix,i,j)
    nc = nblocks(m)
    @boundscheck begin
        checkbounds(Base.OneTo(nc),i)
        checkbounds(Base.OneTo(nc),j)
    end

    @inbounds begin
        return AssocView(m,i,j)
    end
end

@inline function Base.size(m::AssocView)
    _,_,a,b = idx_to_ijab(m.ijab)
    return (a,b)
end

function Base.summary(io::IO,m::AssocView{T}) where T
    s1,s2 = size(m)
    print(io,s1,"×",s2)
    print(io," ")
    i,j,a,b = idx_to_ijab(m.ijab)
    if i == j
        print(io,low_color("(symmetric) "))
    end

    if i > j
        print(io,low_color("(transposed) "))
    end
    print(io,typeof(m))
end

function Base.print_array(io::IO, m::AssocView{T}) where T
    s = size(m)
    M = Matrix{T}(undef,s)
    M .= 0
    for i in 1:s[1]
        for j in 1:s[2]
            M[i,j] = m[i,j]
        end
    end
    Base.print_array(io,M)
end

Base.eltype(m::AssocView{T}) where T = T

#returns the absolute index. that is. it is directly indexable by the parent array
@inline function validindex(m::AssocView{TT},k1::Int,k2::Int) where TT
    ijab = m.ijab
    _i,_j,_a,_b = idx_to_ijab(ijab) #unwrap the ijab, where ij is the specified ij, and ab is the size of the block ij
    (iszero(_a) || iszero(_b)) && return zero(eltype(ijab)) #zero--sized assoc view has no indices
    return validindex(m.m.indices,(_i,_j,k1,k2))
end

@inline validindex(m::Compressed4DMatrix{TT},ijab::NTuple{4,Int}) where TT = validindex(m.indices,ijab)

@inline function validindex(idxs::AbstractVector{T},ijab::NTuple{4,Int}) where T <: Integer
    can = canonical_index(ijab)
    k  = ijab_to_idx(canonical_index(ijab))
    len = length(idxs)
    w = searchsortedfirst(idxs,k,T(1),T(len),Base.Order.ForwardOrdering()) #search unique compressed index
    w > len && return T(0)
    in_idxs = @inbounds(idxs[w]) == k
    return ifelse(in_idxs,T(w),T(0))
end

@inline function Base.getindex(m::AssocView{T},i::Int,j::Int) where T
    idx = validindex(m,i,j)
    iszero(idx) && return _zero(T)
    @inbounds begin
        return m.m.values[idx]
    end
end

@inline function Base.getindex(m::Compressed4DMatrix{T},ijab::NTuple{4,Int}) where T
    idx = validindex(m,ijab)
    vals = m.values
    iszero(idx) && return _zero(T)
    @boundscheck begin
        checkbounds(m.indices,idx)
    end
    @inbounds begin
        return vals[idx]
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

#returns the absolute index. that is. it is directly indexable by the parent array
@inline function validindex_forced!(m::AssocView{TT},k1::Int,k2::Int) where TT
    ijab = m.ijab
    _i,_j,_a,_b = idx_to_ijab(ijab) #unwrap the ijab, where ij is the specified ij, and ab is the size of the block ij
    (iszero(_a) || iszero(_b)) && return zero(eltype(ijab)) #zero--sized assoc view has no indices
    return validindex_forced!(m.m,(_i,_j,k1,k2))
end

@inline function validindex_forced!(m::Compressed4DMatrix{TT},ijab::NTuple{4,Int}) where TT
    idxs = m.indices
    can = canonical_index(ijab)
    k  = ijab_to_idx(canonical_index(ijab))
    len = length(idxs)
    T = eltype(idxs)
    w = searchsortedfirst(idxs,k,T(1),T(len),Base.Order.ForwardOrdering()) #search unique compressed index
    if w > len
        resize!(idxs,len+1)
        resize!(m.values,len+1)
        idxs[end] = k
    end
    in_idxs = @inbounds(idxs[w]) == k
    if !in_idxs
        insert_at!(idxs,w,k)
        insert_at!(m.values,w,_zero(TT))
    end
    return T(w)
end

@inline function Base.setindex!(m::AssocView{T},value,a::Int,b::Int,::Val{true}) where T
    idx = validindex_forced!(m,a,b)
    vals = m.m.values
    @boundscheck begin
        iszero(idx) && throw(BoundsError(m.m, idx))
        checkbounds(m,(a,b))
    end
    @inbounds begin
        vals[idx] = value
    end
end

@inline function Base.setindex!(m::Compressed4DMatrix{T},value,ijab::NTuple{4,Int},::Val{true}) where T
    idx = validindex_forced!(m,ijab)
    vals = m.values
    @boundscheck begin
        iszero(idx) && throw(BoundsError(m, idx))
        i,j,a,b = ijab
        mij = AssocView(m,i,j)
        checkbounds(mij,(a,b))
    end
    @inbounds begin
        vals[idx] = value
    end
end

site_indices(m::Compressed4DMatrix, w::Int) = site_indices(m,idx_to_ijab(m, w))

@inline function site_indices(m::Compressed4DMatrix,ijab::NTuple{4,Int})
    i,j,a,b = ijab
    ia = m.site_offsets[i] + a - 1
    jb = m.site_offsets[j] + b - 1
    return ia, jb
end

function Base.Matrix(m::Compressed4DMatrix{T}) where T
    N = matsize(m)   # total number of sites
    A = fill(_zero(T), (N, N))
    vals = m.values
    for (w, val) in enumerate(vals)
        @inbounds begin
            ia, jb = site_indices(m, w)
            A[ia, jb] = val
            A[jb, ia] = val # symmetry
        end
    end
    return A
end

#compatibility functions

function infer_site_offsets(ij::AbstractVector,ab::AbstractVector,nc = 0)
    if iszero(nc)
        _nc = maximum(maximum,ij)
    else
        _nc = nc
    end
    bsizes = zeros(Int,_nc)
    for k in 1:length(ij)
        i,j = ij[k]
        a,b = ab[k]
        bsizes[i] = max(bsizes[i],a)
        bsizes[j] = max(bsizes[j],b)
    end
    return offsets_from_bsizes!(bsizes)
end

function infer_site_offsets(ijab,nc = 0)
    if iszero(nc)
        _nc = maximum(k -> max(k[1],k[2]),ijab)
    else
        _nc = nc
    end
    bsizes = zeros(Int,_nc)
    for k in 1:length(ijab)
        i,j,a,b = ijab[k]
        bsizes[i] = max(bsizes[i],a)
        bsizes[j] = max(bsizes[j],b)
    end
    return offsets_from_bsizes!(bsizes)
end


function Compressed4DMatrix(vals::AbstractVector{T},ijab::AbstractVector{NTuple{4,Int}},offs::AbstractVector{Int} = infer_site_offsets(ijab)) where {T}
    m = c4d_from_site_offsets(T,offs)
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

function Compressed4DMatrix(vals::AbstractVector{T},ij::Vector{NTuple{2,Int}},ab::Vector{NTuple{2,Int}},offs::AbstractVector{Int} = infer_site_offsets(ij,ab)) where T
    m = c4d_from_site_offsets(T,offs)
    for idx in 1:length(ij)
        _ij = ij[idx]
        _ab = ab[idx]
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
using .Compressed4DMatrices: AssocView, indices, dropzeros!, validindex, idx_to_ijab, ijab_to_idx, assoc_is_pure, assoc_self_assoc, matsize
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
        if ia != jb   # avoid double‑counting diagonal self‑term
            Y[jb] += v * X[ia]
        end
    end
    return Y
end
raw_values(x::Compressed4DMatrix) = x.values

raw_values(x::PackedVofV) = x.v

_param_from_values(x::X,m::Compressed4DMatrix) where X = Compressed4DMatrix(x,m.indices,m.site_offsets)

_param_from_values(x::X,m::P) where {X,P <: PackedVofV} = PackedVofV(m.p,x)

Solvers.primalval(x::T) where T <: Compressed4DMatrix = param_from_values(Solvers.primalval(x.values),x)
Solvers.primalval_eager(x::T) where T <: Compressed4DMatrix = param_from_values(Solvers.primalval_eager(x.values),x)

Solvers.primalval(x::T) where T <: PackedVofV = param_from_values(Solvers.primalval(x.v),x)
Solvers.primalval_eager(x::T) where T <: PackedVofV = param_from_values(Solvers.primalval_eager(x.v),x)

paramtype(::Type{M}) where M <: Compressed4DMatrix{T} where T = T
paramtype(::Type{M}) where M <: PackedVectorsOfVectors.PackedVectorOfVectors{K,V} where {K,V <: AbstractVector{T}} where T = T

function pack_vectors(x::AbstractVector{<:AbstractVector})
    return PackedVectorsOfVectors.pack(x)
end

function packed_zeros!(T,bsizes)
    values = zeros(T,sum(bsizes))
    offsets = Compressed4DMatrices.offsets_from_bsizes!(bsizes)
    return PackedVofV(offsets,values)
end