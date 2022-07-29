import Optim

"""
    isstable(model,V,T,z)::Bool

Performs stability tests for a (V,T,z) pair, and warn if any tests fail. returns `true/false`.

Checks:
 - mechanical stability: isothermal compressibility is not negative.
 - diffusive stability: all eigenvalues of `∂²A/∂n²` are positive.
 
"""
function isstable(model, V, T, z)
    stable = true
    if !mechanical_stability(model, V, T, z)
        @warn "StabilityWarning: Phase is mechanically unstable"
        stable = false
    end
    if !diffusive_stability(model, V, T, z)
        @warn "StabilityWarning: Phase is diffusively unstable"
        stable = false
    end
    if !chemical_stability(model, V, T, z)
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
function mechanical_stability(model, V, T, z)
    return VT_isothermal_compressibility(model, V, T, z) >= 0
end

"""
    diffusive_stability(model,V,T,z)::Bool

Performs a diffusive stability for a (V,T,z) pair, returns `true/false`.
Checks if all eigenvalues of `∂²A/∂n²` are positive.
"""
function diffusive_stability(model, V, T, z)
    isone(length(model)) && return true
    A(x) = eos(model, V, T, x)
    Hf = ForwardDiff.hessian(A, z)
    λ = eigmin(Hf) # calculating just minimum eigenvalue more efficient than calculating all & finding min
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
function gibbs_duhem(model, V, T, z=SA[1.0])
    μ = dot(z, Clapeyron.VT_chemical_potential(model, V, T, z))
    g = VT_gibbs_free_energy(model, V, T, z)
    return abs(μ - g), μ, g
end

"""
    chemical_stability_check(model,p,T,z)::Bool

Performs a chemical stability check using the tangent plane distance criterion, starting with the wilson correlation for K-values.
"""
function chemical_stability(model, V, T, z)
    # in case of only pure components.
    if isone(length(z))
        return pure_chemical_instability(model, V / sum(z), T)
    end

    # For mixtures
    p = pressure(model, V, T, z)
    tm_min_vec, _ = chemical_stability_analysis(model, p, T, z; converge_min=false)

    if minimum(tm_min_vec) .< 0.0
        stable = false
    else
        stable = true
    end

    return stable
end

function chemical_stability_analysis(model, p, T, z; converge_min=true, abstol=1e-3)
    Kʷ = Clapeyron.wilson_k_values(model, p, T)
    z = z ./ sum(z)
    # Generate initial guesses
    w_vap = normalize(z ./ Kʷ, 1) # vapour-like root
    w_liq = normalize(Kʷ .* z, 1) # liquid-like root
    w0_vec = [w_liq, w_vap]

    # Objective function - Unconstrained formulation in mole numbers
    φ(x) = fugacity_coefficient(model, p, T, x)
    d(x) = log.(x) .+ log.(φ(x))
    d_z = d(z)

    tm(W) = 1.0 + sum(W .* (d(W) .- d_z .- 1.0))
    f(W) = tm(exp10.(W))

    if converge_min
        f_callback = (_) -> false
    else
        f_callback = (x) -> f(x) < -eps()
    end

    options = Optim.Options(callback=f_callback, g_tol=abstol)
    sol(w0) = Optim.optimize(f, log10.(w0), Optim.NewtonTrustRegion(), options; autodiff=:forward)

    sol_vec = sol.(w0_vec)
    tm_min_vec = [s.minimum for s in sol_vec]
    tm_xmin_vec = normalize.([exp10.(s.minimizer) for s in sol_vec], 1)
    return tm_min_vec, tm_xmin_vec
end

function pure_chemical_instability(model, V, T)
    Tc, _, _ = crit_pure(model)
    T >= Tc && return true
    psat, vl, vv = saturation_pressure(model, T)
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

export isstable
export mechanical_stability, diffusive_stability, chemical_stability
export gibbs_duhem
