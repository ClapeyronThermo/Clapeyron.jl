# Todo: define mols for z
# Dispatch original function on Numbers
function z_ustrip(z::V) where V <: AbstractVector{U} where U <: AbstractFloat
    return z
end

function z_ustrip(z::V) where V <: AbstractVector{U} where U <: Real
    return float.(z)
end

function z_ustrip(z::V) where V <: AbstractVector{U} where U <: Unitful.Amount
    return float.(ustrip.(u"mol", z))
end

function volume(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"m^3")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, volume(model, p, T, z; phase=phase)*u"m^3")
end

function saturation_pressure(model::SAFT, T::Unitful.Temperature; output=[u"Pa", u"m^3", u"m^3"])
    T = float(ustrip(u"K", T))
    (P_sat, v_l, v_v) = saturation_pressure(model, T)
    return (P_sat*u"Pa", v_l*u"m^3", v_v*u"m^3")
end

function enthalpy_vap(model::SAFT, T::Unitful.Temperature)
    T = float(ustrip(u"K", T))
    return uconvert(output, enthalpy_vap(model, T)*u"J/mol")
end

function pressure(model::SAFT, v::Unitful.Volume, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"Pa")
    z = z_ustrip(z)
    v = float(ustrip(u"m^3", v)); T = float(ustrip(u"K", T))
    return uconvert(output, pressure(model, z, v, T; phase=phase)*u"Pa")
end

function entropy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, entropy(model, p, T, z; phase=phase)*u"J")
end

function chemical_potential(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, chemical_potential(model, p, T, z; phase=phase)*u"J")
end

function internal_energy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, internal_energy(model, p, T, z; phase=phase)*u"J")
end

function enthalpy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, enthalpy(model, p, T, z; phase=phase)*u"J")
end

function gibbs_free_energy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, gibbs_free_energy(model, p, T, z; phase=phase)*u"J")
end

function helmholtz_free_energy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, helmholtz_free_energy(model, p, T, z; phase=phase)*u"J")
end

function isochoric_heat_capacity(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J/K")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, isochoric_heat_capacity(model, p, T, z; phase=phase)*u"J/K")
end

function isobaric_heat_capacity(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J/K")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, isobaric_heat_capacity(model, p, T, z; phase=phase)*u"J/K")
end

function isothermal_compressibility(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"1/Pa")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, isothermal_compressibility(model, p, T, z; phase=phase)*u"1/Pa")
end

function isentropic_compressibility(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"1/Pa")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, isentropic_compressibility(model, p, T, z; phase=phase)*u"1/Pa")
end

function speed_of_sound(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"m/s")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, speed_of_sound(model, p, T, z; phase=phase)*u"m/s")
end

function isobaric_expansivity(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"1/K")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, isobaric_expansivity(model, p, T, z; phase=phase)*u"1/K")
end

function joule_thomson_coefficient(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"K/Pa")
    z = z_ustrip(z)
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, joule_thomson_coefficient(model, p, T, z; phase=phase)*u"K/Pa")
end
