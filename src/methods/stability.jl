"""
    isstable(model,p,T,z)::Bool

Performs stability tests for a (p,T,z) pair, and warn if any tests fail. Returns `true/false`.

Checks, in order of complexity:
 - mechanical stability: isothermal compressibility is not negative.
 - diffusive stability: all eigenvalues of `∂²A/∂n²` are positive.
 - chemical stability: there isn't any other combinations of compositions at p(V,T),T that are more stable than the input composition.

For checking (V,T,z) pairs, use `Clapeyron.VT_isstable(model,V,T,z)` instead.
"""
function isstable(model,p,T,z = SA[1.0])
    V = volume(model,p,T,z)
    return VT_isstable(model,V,T,z,false)
end

function VT_isstable(model,V,T,z = SA[1.0],check_v = true)
    stable = true
    VT_mechanical_stability(model,V,T,z) || return false
    VT_diffusive_stability(model,V,T,z) || return false
    VT_chemical_stability(model,V,T,z,check_v) || return false
    return true
end

"""
    VT_mechanical_stability(model,V,T,z = SA[1.0])::Bool

Performs a mechanical stability for a (V,T,z) pair, returns `true/false`.
Checks if isothermal compressibility is not negative.


!!! note
    This function does not have a `p`,`T` counterpart, because if we calculate the volume via `volume(model,p,T,z)`, it will be, by definition, a mechanically stable phase.
"""
function VT_mechanical_stability(model,V,T,z = SA[1.0])
    return VT_isothermal_compressibility(model,V,T,z) >= 0
end

"""
    VT_diffusive_stability(model,V,T,z)::Bool

Performs a diffusive stability for a (V,T,z) pair, returns `true/false`.
Checks if all eigenvalues of `∂²A/∂n²` are positive.
Returns `false` if the eos calculation failed. Rhis normally occurs when evaluating on densities lower than the maximum density (given by `Clapeyron.lb_volume(model,T,z)`).
"""
function VT_diffusive_stability(model,V,T,z = SA[1.0])
    ρᵢ = similar(z,Base.promote_eltype(V,z))
    ρᵢ .= z ./ V
    HΨ = Ψ_hessian(model,T,ρᵢ)
    if any(!isfinite,HΨ)
        return false
    end
    if length(model) == 1
        return HΨ[1,1] > 0
    end
    λ = eigmin(Hermitian(HΨ)) # calculating just minimum eigenvalue more efficient than calculating all & finding min
    return λ > 0
end

function VT_diffusive_eigvalue(model,V,T,z = SA[1.0])
    ρᵢ = similar(z,Base.promote_eltype(V,z))
    ρᵢ .= z ./ V
    HΨ = Ψ_hessian(model,T,ρᵢ)
    if any(!isfinite,HΨ)
        return zero(eltype(HΨ))/zero(eltype(HΨ))
    end
    if length(model) == 1
        return HΨ[1,1]
    end
    return eigmin(Hermitian(HΨ)) # calculating just minimum eigenvalue more efficient than calculating all & finding min
end

"""
    diffusive_stability(model,p,T,z = SA[1.0],phase = :unknown,threaded = true,vol0 = nothing)

Performs a diffusive stability for a (V,T,z) pair, returns `true/false`.

    Checks if all eigenvalues of `∂²A/∂n²` are positive.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

"""
function diffusive_stability(model,p,T,z = SA[1.0],phase = :unknown,threaded = true,vol0 = nothing)
    V = volume(model,p,T,z;phase,threaded,vol0)
    return VT_diffusive_stability(model,V,T,z)
end

"""
    gibbs_duhem(model,V,T,z=[1.0])

Performs a Gibbs-Duhem check on the input conditions:

```
∑zᵢμᵢ - G ≈ 0
```
Where `G` is the total gibbs free energy. It can help diagnose if a user-defined eos is consistent.

Returns |∑zᵢμᵢ - G|, ∑zᵢμᵢ and G at the specified conditions.
"""
function gibbs_duhem(model,V,T,z=SA[1.0])
    μ = dot(z,Clapeyron.VT_chemical_potential(model,V,T,z))
    g = VT_gibbs_free_energy(model,V,T,z)
    return abs(μ-g),μ,g
end

"""
    ideal_consistency(model,V,T,z=[1.0])

Performs a ideal model consistency check:
```
∂a₀∂V + 1/V ≈ 0
```
Where `∂a₀∂V` is the derivative of `a_ideal` respect to `V`. It can help diagnose if a user-defined ideal model is consistent.

Return |∂a₀∂V + 1/V| at the specified conditions.

If the model is not an `IdealModel`, then `Clapeyron.idealmodel(model)` will be called to obtain the respective ideal model.
"""

function ideal_consistency(model,V,T,z =SA[1.0])
    id = idealmodel(model)
    if id === nothing || id === model
        f(∂V) = a_ideal(model,∂V,T,z)
        ∂f0∂V = Solvers.derivative(f,V)
        n = sum(z)
        return abs(∂f0∂V + 1/V)
    else
        return ideal_consistency(id,V,T,z)
    end
end

"""
    VT_chemical_stability(model,V,T,z)::Bool

Performs a chemical stability check using the tangent plane distance criterion, using the [tpd](@ref) function
"""
function VT_chemical_stability(model::EoSModel,V,T,z,check_vol = true)
    if isone(length(z))
        return pure_chemical_instability(model,V/sum(z),T)
    end
    p = pressure(model,V,T,z)
    if check_vol
        Vx = volume(model,p,T,z)
        #we check first if the phase itself is stable, maybe there is another phase
        #with the same composition, but with a different volume, that is more stable.
        V ≈ Vx || return false
    end
    #the input volume corresponds to the most stable phase at that compsition.
    #proceed to tpd
    return chemical_stability(model,p,T,z)
end

function chemical_stability(model,p,T,z)
    length(model) == 1 && return true #there aren't other combinations of composition.
    comps,wi,_,_ = tpd(model,p,T,z,break_first = true)
    return iszero(length(comps))
end

function pure_chemical_instability(model,V,T)
    Tc,_,_ = crit_pure(model)
    T >= Tc && return true
    psat,vl,vv = saturation_pressure(model,T)
    if isnan(psat)
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
    μ(w) = VT_chemical_potential(model,V,T,w)
    tdp = sum(w.*(μ(w) .- μ(z)))#./(8.314*T)
end

export isstable
export VT_mechanical_stability, VT_diffusive_stability,VT_chemical_stability
export diffusive_stability, chemical_stability
export gibbs_duhem
