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
    if !chemical_stability(model,V,T,z)
        @warn "StabilityWarning: Phase is chemically unstable"
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

"""
    chemical_stability(model,V,T,z)::Bool
Performs a chemical stability check using the 
"""
function chemical_stability(model::EoSModel,p,T,z)
    # Generate vapourlike and liquidlike initial guesses
    # Currently using Wilson correlation
    Pc = model.params.Pc.values
    Tc = model.params.Tc.values
    ω = model.alpha.params.acentricfactor.values

    Kʷ = @. Pc/p*exp(5.373*(1+ω)*(1-Tc/T))
    z = z./sum(z)
    w_vap = Kʷ.*z
    w_liq = z./Kʷ

    tdp_func(w) = Optim.minimum(optimize(w -> tangent_plane_distance(model,p,T,z,w), z))
    tdp = tdp_func(w_vap), tdp_func(w_liq)
    if any(tdp > 0)
        return true
    else
        return false
    end
end

"""
    tangent_plane_distance(model,V,T,z)::Float
Calculates the tangent plane distance for a tangent plane stability test
Uses unconstrained minimisation from NLSolvers.jl
"""
function tangent_plane_distance(model,p,T,z,w)
    w = w./sum(w)
    V = volume(model, p, T, w)

    μ(w) = Clapeyron.VT_chemical_potential(model,V,T,w)

    tdp = sum(w.*(μ(w) .- μ(z)))./(8.314*T)
end
