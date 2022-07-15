"""
    isstable(model,V,T,z)::Bool

Performs stability tests for a (V,T,z) pair, and warn if any tests fail. returns `true/false`.

Checks:
 - mechanical stability: isothermal compressibility is not negative.
 - diffusive stability: all eigenvalues of `∂²A/∂n²` are positive.
 
"""
function isstable(model, V, T, z, phase=:stable)
    if phase != :stable
        return true
    end
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
    ((tm_min, λ_min), _) = chemical_stability_analysis(model, p, T, z)

    stable = true
    if λ_min < 0 # Unstable
        stable = false
    elseif tm_min < 0 # Metastable
        stable = false #! How should a metastable state be reported?
    end
    return stable
end

# Michelsen p,T stability analysis
function chemical_stability_analysis(model::EoSModel, p, T, z)
    # Generate vapourlike and liquidlike initial guesses
    # Currently using Wilson correlation

    Kʷ = wilson_k_values(model, p, T)
    z = z ./ sum(z)
    w_liq = z ./ Kʷ
    w_vap = Kʷ .* z

    # Objective function - Unconstrained formulation in mole numbers
    φ(x) = fugacity_coefficient(model, p, T, x)
    d(x) = log.(x) .+ log.(φ(x))
    tm(W, Z) = 1.0 + sum(W .* (d(W) .- d(Z) .- 1))

    # This currently uses Newton, BFGS performs better in my tests with Optim
    # Also would ideally be solved with far lower tolerance - stability analysis doesn't have to converge to machine precision
    tm_func(W) = Solvers.optimize(W -> tm(W .^ 2, z), sqrt.(W))
    res_vec = [tm_func(w_liq), tm_func(w_vap)]
    (tm_min, idx) = findmin(map(x -> x.info.minimum, res_vec))
    tm_xmin = map(x -> sqrt.(abs.(x.info.solution)), res_vec)[idx]
    phase = (idx == 1) ? "liq" : "vap"
    println("stable phase = $phase")
    @show round.(tm_xmin, digits=2)
    return map(x -> x.info.minimum, res_vec)
    # We only need to consider the Hessian if wanting to distinguish between Metastable and unstable states
    # H = ForwardDiff.hessian(W -> tm(W, z), tm_xmin)
    # H_scaled = zeros(size(H))
    # for i in 1:size(H)[1]
    #     for j in 1:size(H)[2]
    #         H_scaled[i, j] = sqrt(tm_xmin[i] * tm_xmin[j]) * H[i, j]
    #     end
    # end
    # λ_min = eigmin(H_scaled)

    # Return compositions for initial guesses in flash algorithms
    # return ((tm_min, 1.0), normalize(tm_xmin, 1))
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
