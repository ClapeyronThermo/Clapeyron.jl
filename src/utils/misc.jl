"""
    gibbs_duhem(model,V,T,z=[1.0])

performs a Gibbs-Duhem check on the input conditions:

```
∑zᵢμᵢ - G ≈ 0
```
Where `G` is the total gibbs energy. it can help diagnose if a user-defined eos is consistent.

return |∑zᵢμᵢ - G|, ∑zᵢμᵢ and G at the specified conditions.
"""
function gibbs_duhem(model,V,T,z=SA[1.0])
    μ = dot(z,Clapeyron.VT_chemical_potential(model,V,T,z))
    g = VT_gibbs_free_energy(model,V,T,z)
    return abs(μ-g),μ,g
end