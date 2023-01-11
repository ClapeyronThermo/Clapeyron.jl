"""
    isstable(model,V,T,z)::Bool

Performs stability tests for a (V,T,z) pair, and warn if any tests fail. returns `true/false`.

Checks:
 - mechanical stability: isothermal compressibility is not negative.
 - diffusive stability: all eigenvalues of `∂²A/∂n²` are positive.
 
"""
function isstable(model,V,T,z)
    stable = true
    if !VT_mechanical_stability(model,V,T,z)
        @warn "StabilityWarning: Phase is mechanically unstable"
        stable = false
    end
    if !VT_diffusive_stability(model,V,T,z)
        @warn "StabilityWarning: Phase is diffusively unstable"
        stable = false
    end
    if !VT_chemical_stability(model,V,T,z)
        @warn "StabilityWarning: Phase is chemically unstable"
        stable = false
    end
    return stable
end

"""
    VT_mechanical_stability(model,V,T,z)::Bool

Performs a mechanical stability for a (V,T,z) pair, returns `true/false`.
Checks if isothermal compressibility is not negative. 
"""
function VT_mechanical_stability(model,V,T,z)
    return VT_isothermal_compressibility(model,V,T,z) >= 0
end

"""
    VT_diffusive_stability(model,V,T,z)::Bool

Performs a diffusive stability for a (V,T,z) pair, returns `true/false`.
Checks if all eigenvalues of `∂²A/∂n²` are positive.
"""
function VT_diffusive_stability(model,V,T,z)
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

"""
    VT_chemical_stability(model,V,T,z)::Bool

Performs a chemical stability check using the tangent plane distance criterion, starting with the wilson correlation for K-values.
"""
function VT_chemical_stability(model::EoSModel,V,T,z)
    # Generate vapourlike and liquidlike initial guesses
    # Currently using Wilson correlation

    #in case of only pure components.
    if isone(length(z))
        return pure_chemical_instability(model,V/sum(z),T) 
    end
    p = pressure(model,V,T,z)
    Kʷ = wilson_k_values(model,p,T)
    z = z./sum(z)
    w_vap = Kʷ.*z
    w_liq = z./Kʷ
   
    tdp_func(w,phase) = Solvers.optimize(w -> tangent_plane_distance(model,p,T,z,phase,w), w) |> Solvers.x_minimum
    tdp = (tdp_func(w_vap,:v), tdp_func(w_liq,:l))
    if minimum(tdp) >= 0
        return true
    else
        return false
    end
end

function pure_chemical_instability(model,V,T)
    Tc,Pc,Vc = crit_pure(model)
    T >= Tc && return true
    psat,vl,vv = saturation_pressure(model,T)
    if isnan(psat)
        @error "could not determine chemical instability. saturation solver failed."
        return false
    end
    if (vl < V < vv)
        return false
    else
        return true
    end
end

"""
    tangent_plane_distance(model,V,T,z,phase::Symbol,w)::Float
Calculates the tangent plane distance for a tangent plane stability test
"""
function tangent_plane_distance(model,p,T,z,phase,w)
    w = w./sum(w)
    V = volume(model, p, T, w;phase=phase)
    μ(w) = Clapeyron.VT_chemical_potential(model,V,T,w)
    tdp = sum(w.*(μ(w) .- μ(z)))#./(8.314*T)
end

export isstable
export VT_mechanical_stability, VT_diffusive_stability,VT_chemical_stability
export gibbs_duhem
