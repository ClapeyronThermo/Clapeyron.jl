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
        #@warn "StabilityWarning: Phase is mechanically unstable"
        stable = false
    end
    if !VT_diffusive_stability(model,V,T,z)
        #@warn "StabilityWarning: Phase is diffusively unstable"
        stable = false
    end
    if !chemical_stability(model,V,T,z)
        #@warn "StabilityWarning: Phase is chemically unstable"
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
function gibbs_duhem(model,V,T,z=SA[1.0])
    μ = dot(z,Clapeyron.VT_chemical_potential(model,V,T,z))
    g = VT_gibbs_free_energy(model,V,T,z)
    return abs(μ-g),μ,g
end

"""
    ideal_consistency(model,V,T,z=[1.0])

performs a ideal model consistency check:
```
∂a₀∂V + 1/V ≈ 0
```
Where `∂a₀∂V` is the derivative of `a_ideal` respect to `V`. it can help diagnose if a user-defined ideal model is consistent.

Return |∂a₀∂V + 1/V| at the specified conditions.

If the model is not an `IdealModel`, then `Clapeyron.idealmodel(model)` will be called to obtain the respective ideal model.
"""

function ideal_consistency(model,V,T,z =SA[1.0])
    id = idealmodel(model)
    if id === nothing
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
    Kʷ = tp_flash_K0(model,p,T)
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
    _, tm_min_vec = chemical_stability_analysis(model, p, T, z; converge_min=false)

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
    _φ(x) = fugacity_coefficient(model, p, T, x)
    _d(x) = log.(x) .+ log.(_φ(x))
    d_z = _d(z)

    function f(logW)
        W = exp.(logW)
        φ = fugacity_coefficient(model, p, T, W)
        φ .= log.(W) .+ log.(φ)
        d_W = φ
        tm = 1.0 + @sum(W[i]*(d_W[i] - d_z[i] - 1.0))
        return tm
    end
    #=
    if converge_min
        f_callback = (_) -> false
    else
        f_callback = (x) -> f(x) < -eps()
    end =#

    #options = Optim.Options(callback=f_callback, g_tol=abstol)
    sol(w0) = Solvers.optimize(f, log.(w0))#, options)

    sol_vec = sol.(w0_vec)
    tm_min_vec = [s.minimum for s in sol_vec]
    tm_xmin_vec = normalize.([exp.(s.minimizer) for s in sol_vec], 1)
    return tm_xmin_vec, tm_min_vec
end

function pure_chemical_instability(model,V,T)
    Tc,_,_ = crit_pure(model)
    T >= Tc && return true
    psat,vl,vv = saturation_pressure(model,T)
    if isnan(psat)
        #@error "could not determine chemical instability. saturation solver failed."
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
