
function entropy(model::EoSModel, p, T,   z=SA[1.]; phase = :unknown)
    v      = volume(model, p, T, z; phase=phase)
    return vt_entropy(model,v,T,z)
end

function chemical_potential(model::EoSModel, p, T, z= SA[1.]; phase = :unknown)
    v      = volume(model, p, T, z; phase=phase)
    return vt_chemical_potential(modle,v,T,z)
end

function internal_energy(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v      = volume(model, p, T, z; phase=phase)
    return vt_internal_energy(model,v,T,z)
end

function enthalpy(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v      = volume(model, p, T, z; phase=phase)
    return vt_enthalpy(model,v,T,z)
end

function gibbs_free_energy(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v      = volume(model, p, T, z; phase=phase)
    return vt_gibbs_free_energy(model,v,T,z)
end

function helmholtz_free_energy(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v      = volume(model, p, T, z; phase=phase)
    return vt_helmholtz_free_energy(model,v,T,z)
end

function isochoric_heat_capacity(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v       = volume(model, p, T, z; phase=phase)
    return vt_isochoric_heat_capacity(model,v,T,z)
end

function isobaric_heat_capacity(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v       = volume(model, p, T, z; phase=phase)
    return vt_isobaric_heat_capacity(model,v,T,z)
end

function isothermal_compressibility(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v       = volume(model, p, T, z; phase=phase)
    return vt_isothermal_compressibility(model,v,T,z)
end

function isentropic_compressibility(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v       = volume(model, p, T, z; phase=phase)
    return vt_isentropic_expansivity(model,v,T,z)
end

function speed_of_sound(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v       = volume(model, p, T, z; phase=phase)
    return vt_speed_of_sound(model,v,T,z)
end

function isobaric_expansivity(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v = volume(model, p, T, z; phase=phase)
    return vt_isobaric_expansivity(model,v,T,z)
end

function joule_thomson_coefficient(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v  = volume(model, p, T, z; phase=phase)
    return vt_joule_thomson_coefficient(model,v,T,z)
end


function compressibility_factor(model::EoSModel, p, T,  z=SA[1.]; phase = :unknown)
    v  = volume(model, p, T, z; phase=phase)
    return p*v/(R̄*T)
end

function inversion_temperature(model::EoSModel,p,z=SA[1.0])
    T0 = 6.75*T_scale(model,z)
    μⱼₜ(T) = joule_thomson_coefficient(model,p,T,z)
    return Roots.find_zero(μⱼₜ,T0)
end

export entropy, chemical_potential, internal_energy, enthalpy, gibbs_free_energy
export helmholtz_free_energy, isochoric_heat_capacity, isobaric_heat_capacity
export isothermal_compressibility, isentropic_compressibility, speed_of_sound
export isobaric_expansivity, joule_thomson_coefficient, compressibility_factor, inversion_temperature
