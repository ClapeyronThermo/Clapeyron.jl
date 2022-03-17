"""
    isstable(model,V,T,z)::Bool

Performs stability tests for a (V,T,z) pair, and warn if any tests fail. returns `true/false`.

Checks:
 - mechanical stability: isothermal compressibility is not negative.
 - diffusive stability: all eigenvalues of `∂²A/∂n²` are positive.
 
"""
function isstable(model,V,T,z,phase=:stable)
    if phase != :stable
        return true
    end
    stable = true
    if !mechanical_stability(model,V,T,z)
        @warn "StabilityWarning: Phase is mechanically unstable"
        stable = false
    end
    if !diffusive_stability(model,V,T,z)
        @warn "StabilityWarning: Phase is diffusively unstable"
        stable = false
    end
    return stable
end

"""
    mechanical_stability(model,V,T,z)::Bool

Performs a mechanical stability for a (V,T,z) pair, returns `true/false`.
Checks if isothermal compressibility is not negative. 
"""
function mechanical_stability(model,V,T,z)
    return VT_isothermal_compressibility(model,V,T,z) >= 0
end

"""
    diffusive_stability(model,V,T,z)::Bool

Performs a diffusive stability for a (V,T,z) pair, returns `true/false`.
Checks if all eigenvalues of `∂²A/∂n²` are positive.
"""
function diffusive_stability(model,V,T,z)
    isone(length(model)) && return true
    A(x) = eos(model,V,T,x)
    Hf = ForwardDiff.hessian(A,z)
    (Λ,U)=eigen(Hf)
    λ = minimum(Λ)
    return λ > 0
end

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

export isstable, mechanical_stability, diffusive_stability, gibbs_duhem


