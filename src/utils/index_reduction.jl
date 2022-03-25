"""
    index_reduction(model::EoSModel,z,zmin = sum(z)*4*eps(eltype(z)))

Removes any component with composition `z[i] < zmin`. returns a reduced model `model_r` and a vector of indices `idx_r`, such as:

```julia
model_r,idx_r = index_reduction(model,z)
eos(model,V,T,z) ≈ eos(model_r,V,T,z[idx_r])
```
if the model does not have empty compositions, it will just return the input model
"""
function index_reduction(model,z,zmin = sum(z)*4*eps(eltype(z)))
    length(model) == 1 && (return model,[1])
    idx_r = findall(>(zmin),z)
    if length(idx_r) == length(z)
        model_r = model
    else
        model_r = split_model(model,[idx_r]) |> only
    end
    return model_r,idx_r
end

export index_reduction