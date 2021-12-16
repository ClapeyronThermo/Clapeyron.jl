function entropy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_entropy(model,V,T,z)
end

function entropy_res(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_entropy_res(model,V,T,z)
end

function chemical_potential(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_chemical_potential(model,V,T,z)
end

function chemical_potential_res(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_chemical_potential_res(model,V,T,z)
end

function internal_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_internal_energy(model,V,T,z)
end

function enthalpy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_enthalpy(model,V,T,z)
end

function gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_gibbs_free_energy(model,V,T,z)
end

function helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_helmholtz_free_energy(model,V,T,z)
end

function isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isochoric_heat_capacity(model,V,T,z)
end

function isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isobaric_heat_capacity(model,V,T,z)
end

function isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isothermal_compressibility(model,V,T,z)
end

function isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isentropic_compressibility(model,V,T,z)
end

function speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_speed_of_sound(model,V,T,z)
end

function isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isobaric_expansivity(model,V,T,z)
end

function joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_joule_thomson_coefficient(model,V,T,z)
end

function fugacity_coefficient(model::EoSModel,p,T,z=SA[1.]; phase = :unknown, threaded=true)
    V = volume(model,p,T,z;phase,threaded)
    μ_res = VT_chemical_potential_res(model,V,T,z)
    φ = μ_res
    Z = p*V/R̄/T/sum(z)
    for i ∈ @comps
        φ[i] = exp(μ_res[i]/R̄/T)/Z
    end
    return φ
end

function activity_coefficient(model::EoSModel,p,T,z=SA[1.]; phase = :unknown, threaded=true)
    pure   = split_model(model)
    μ_mixt = chemical_potential(model,p,T,z;phase,threaded)
    γ_i = μ_mixt
    for i ∈ @comps
        μ_pure_i = chemical_potential(pure[i],p,T;phase,threaded)[1]
        γ_i[i] = exp((μ_mixt[i]-μ_pure_i)/R̄/T)/z[i]
    end
    return γ_i
end
"""
    compressibility_factor(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)

Calculates the compressibility factor `Z`, defined as:

```julia
Z = p*V(p)/R*T
```
the keywords `phase` and `threaded` are passed to the [volume solver](@ref Clapeyron.volume).
"""
function compressibility_factor(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_compressibility_factor(model, V, T, z)
end

function inversion_temperature(model::EoSModel, p, z=SA[1.0])
    T0 = 6.75*T_scale(model,z)
    μⱼₜ(T) = joule_thomson_coefficient(model,p,T,z)
    return Roots.find_zero(μⱼₜ,T0)
end

function molar_density(model::EoSModel,p,T,z=SA[1.0];phase = :unknown,threaded=true)
     V = volume(model,p,T,z;phase=phase,threaded=threaded)
     return sum(z)/V
end

function mass_density(model::EoSModel,p,T,z=SA[1.0];phase = :unknown,threaded=true)
    V = volume(model,p,T,z;phase=phase,threaded=threaded)
    molar_weight = molecular_weight(model,z)
    return molar_weight/V
end

export entropy, chemical_potential, internal_energy, enthalpy, gibbs_free_energy
export helmholtz_free_energy, isochoric_heat_capacity, isobaric_heat_capacity
export isothermal_compressibility, isentropic_compressibility, speed_of_sound
export isobaric_expansivity, joule_thomson_coefficient, compressibility_factor, inversion_temperature
export mass_density,molar_density, activity_coefficient, fugacity_coefficient, entropy_res