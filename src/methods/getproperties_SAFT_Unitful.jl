# Todo: define mols for z
# Dispatch original function on Numbers

function get_volume(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(Pa, p); T = ustrip(K, T)
    return get_volume(model::SAFT, z, p, T, phase)*u"m^3"
end

function get_sat_pure(model::SAFT, T::Unitful.Temperature)
    T = ustrip(u"K", T)
    (P_sat, v_l, v_v) = get_sat_pure(model, T)
    return (P_sat*u"Pa", v_l*u"m^3", v_v*u"m^3")
end

function get_enthalpy_vap(model::SAFT, T::Unitful.Temperature)
    T = ustrip(u"K", T)
    return get_enthalpy_vap(model, T)*u"J/mol"
end

function get_pressure(model::SAFT, z, v::Unitful.Volume, T::Unitful.Temperature, phase="unknown")
    v = ustrip(u"m^3", v); T = ustrip(u"K", T)
    return get_pressure(model, z, v, T, phase=phase)*u"Pa"
end

function get_entropy(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_entropy(model, z, p, T, phase=phase)*u"J"
end

function get_chemical_potential(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_chemical_potential(model, z, p, T, phase)*u"J"
end

function get_internal_energy(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_internal_energy(model, z, p, T, phase)*u"J"
end

function get_enthalpy(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_enthalpy(model, z, p, T, phase)*u"J"
end

function get_Gibbs_free_energy(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_Gibbs_free_energy(model, z, p, T, phase)*u"J"
end

function get_Helmholtz_free_energy(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_Helmholtz_free_energy(model, z, p, T, phase)*u"J"
end

function get_isochoric_heat_capacity(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_isochoric_heat_capacity(model, z, p, T, phase)*u"J/K"
end

function get_isobaric_heat_capacity(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_isobaric_heat_capacity(model, z, p, T, phase)*u"J/K"
end

function get_isothermal_compressibility(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_isothermal_compressibility(model, z, p, T, phase)*u"1/Pa"
end

function get_isentropic_compressibility(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_isentropic_compressibility(model, z, p, T, phase)*u"1/Pa"
end

function get_speed_of_sound(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_speed_of_sound(model, z, p, T, phase)*u"m/s"
end

function get_isobaric_expansivity(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_isobaric_expansivity(model, z, p, T, phase)*u"1/K"
end

function get_Joule_Thomson_coefficient(model::SAFT, z, p::Unitful.Pressure, T::Unitful.Temperature, phase="unknown")
    p = ustrip(u"Pa", p); T = ustrip(u"K", T)
    return get_Joule_Thomson_coefficient(model, z, p, T, phase)*u"K/Pa"
end
