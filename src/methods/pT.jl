function entropy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_entropy(model,V,T,z)
end

function chemical_potential(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_chemical_potential(model,V,T,z)
end

function internal_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_internal_energy(model,V,T,z)
end

function enthalpy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_enthalpy(model,V,T,z)
end

function gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_gibbs_free_energy(model,V,T,z)
end

function helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_helmholtz_free_energy(model,V,T,z)
end

function isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_isochoric_heat_capacity(model,V,T,z)
end

function isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_isobaric_heat_capacity(model,V,T,z)
end

function isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_isothermal_compressibility(model,V,T,z)
end

function isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_isobaric_expansivity(model,V,T,z)
end

function speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_speed_of_sound(model,V,T,z)
end

function isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_isobaric_expansivity(model,V,T,z)
end

function joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return VT_joule_thomson_coefficient(model,V,T,z)
end

function compressibility_factor(model::EoSModel, p, T, z=SA[1.]; phase = :unknown)
    V = volume(model, p, T, z; phase=phase)
    return p*V/(R̄*T)
end

function inversion_temperature(model::EoSModel, p, z=SA[1.0])
    T0 = 6.75*T_scale(model,z)
    μⱼₜ(T) = joule_thomson_coefficient(model,p,T,z)
    return Roots.find_zero(μⱼₜ,T0)
end

export entropy, chemical_potential, internal_energy, enthalpy, gibbs_free_energy
export helmholtz_free_energy, isochoric_heat_capacity, isobaric_heat_capacity
export isothermal_compressibility, isentropic_compressibility, speed_of_sound
export isobaric_expansivity, joule_thomson_coefficient, compressibility_factor, inversion_temperature
