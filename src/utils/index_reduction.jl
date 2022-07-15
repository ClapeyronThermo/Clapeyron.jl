"""
    index_reduction(model::EoSModel,z,zmin = sum(z)*4*eps(eltype(z)))

Removes any component with composition `z[i] < zmin`. returns a reduced model `model_r` and a vector of indices `idx_r`, such as:

```julia
model_r,idx_r = index_reduction(model,z)
eos(model,V,T,z) ≈ eos(model_r,V,T,z[idx_r])
```
if the model does not have empty compositions, it will just return the input model
"""
function index_reduction(model::EoSModel,z::AbstractVector,zmin = sum(z)*4*eps(eltype(z)))
    length(model) == 1 && (return model,[true])
    idx = z .> zmin
    if all(idx)
        model_r = model
    else
        model_r = split_model(model,[findall(idx)]) |> only
    end
    return model_r,idx
end



function index_expansion(x::Matrix,idr::AbstractVector)
    numspecies = length(idr)
    l1,_ = size(x)
    res = zeros(eltype(x),l1, numspecies)
    res .= 0
    for i in 1:l1
        res[i,idr] .= x[i,:]
    end
    return res
end

function index_expansion(x::AbstractMatrix,idr::AbstractVector)
    numspecies = length(idr)
    l1,_ = size(x)
    res = similar(x,(l1, numspecies))
    res .= 0
    for i in 1:l1
        res[i,idr] .= x[i,:]
    end
    return res
end

"""
    index_expansion(x::Vector,idx::Vector{Bool})

Given an input vector generated from a reduced model and the non zero indices, returns a resized Vector corresponding to the original model.
"""
function index_expansion(x::AbstractVector,idr::AbstractVector)
    numspecies = length(idr)
    res = similar(x, numspecies)
    res .= 0
    res[idr] .= x
    return res
end

export index_reduction