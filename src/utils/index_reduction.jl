"""
    index_reduction(model::EoSModel,z,zmin = sum(z)*4*eps(eltype(z)))
    index_reduction(model::EoSModel,bools <: AbstractVector{Bool})
Removes any component with composition `z[i] < zmin`. returns a reduced model `model_r` and a vector of indices `idx_r`, such as:

```julia
model_r,idx_r = index_reduction(model,z)
eos(model,V,T,z) ≈ eos(model_r,V,T,z[idx_r])
```
if the model does not have empty compositions, it will just return the input model.

The function will error if the reduction results in an empty model.

You can pass an arbitrary boolean vector (`bools`) to perform the reduction.
"""
function index_reduction(model::EoSModel,z::AbstractVector,zmin = sum(abs,z)*4*eps(float(oneunit(eltype(z)))), check = true)
    #skip splitting if possible
    if check
        nc = length(model)
        nz = length(z)
        nz > nc && throw(BoundsError(model,nz))
        nc > nz && throw(BoundsError(z,nc))
        nc == 1 && return model,BitVector((true,))
        nzero = count(x -> (abs(x) <= zmin) ,z)
        nneg = count(x -> x <= -zmin,z)
        nzero == nc && __index_reduction_empty_error(z)
        nneg > 0 && __index_reduction_negative_error(z)
        nzero == 0 && return model,trues(length(model))
        
    end
    idx = BitVector(z[i] > zmin  for i in eachindex(z))
    return index_reduction(model,idx;check = false)
end

#in case someone uses a number instead of a vector
function index_reduction(model::EoSModel,z::Number,zmin = sum(z)*4*eps(eltype(z));check = true)
    nc = length(model)
    nc != 1 && throw(BoundsError(z,nc))
    nx = abs(z) > zmin
    nx == 0 && __index_reduction_empty_error(z)
    return model,BitVector((true,))
end

function index_reduction(model::EoSModel,bools::T;check = true) where T<:AbstractVector{Bool}
    return index_reduction(model,BitVector(bools);check)
end

#to address arbitrary reduction techniques
function index_reduction(model::EoSModel,bools::BitVector;check = true)
    if check
        nc = length(model)
        nz = length(bools)
        nzero = nc - count(bools)
        nz > nc && throw(BoundsError(model,nz))
        nc > nz && throw(BoundsError(bools,nc))
        nzero == nc && __index_reduction_empty_error(bools)
        nc == 1 && return model,BitVector((true,))
        nzero == 0 && return model,trues(length(model))
    end
    return each_split_model(model,bools),bools
end

@noinline function __index_reduction_empty_error(z)
    throw(error("Clapeyron.index_reduction(model,z): inpút composition $z results in an empty model."))
end

@noinline function __index_reduction_negative_error(z)
    throw(error("Clapeyron.index_reduction(model,z): inpút composition $z has significant negative values."))
end

@noinline function __index_reduction_empty_error(z::AbstractVector{Bool})
    throw(error("Clapeyron.index_reduction(model,bools::AbstractVector{Bool}): input reduction vector `bools` results in an empty model."))
end

index_reduction(::Nothing,idr::AbstractVector) = nothing
index_reduction(x::AbstractVector,idr::AbstractVector) = x[idr]

function index_expansion(x::AbstractMatrix,idr::AbstractVector)
    numspecies = length(idr)
    l1,l2 = size(x)
    if l2 == numspecies
        return x
    end
    res = similar(x,(l1, numspecies))
    res .= 0
    for i in 1:l1
        res[i,idr] .= @view(x[i,:])
    end
    return res
end

"""
    index_expansion(x::Vector,idx::Vector{Bool})
    index_expansion(x::Matrix,idx::Vector{Bool})


Given an input vector generated from a reduced model and the non zero indices, returns a Vector corresponding to the original model.
If the sizes of `x` and `idx` are the same, return the original input.
"""
function index_expansion(x::AbstractVector,idr::AbstractVector)
    numspecies = length(idr)
    res = similar(x, numspecies)
    if length(x) == numspecies
        res .= x
        return res
    end
    res .= false
    res[idr] .= x
    return res
end

export index_reduction
