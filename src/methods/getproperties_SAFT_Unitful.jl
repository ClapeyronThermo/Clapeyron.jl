# Todo: define mols for z
# Dispatch original function on Numbers

function get_volume(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"m^3")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_volume(model, p, T, z; phase=phase)*u"m^3")
end

function get_sat_pure(model::SAFT, T::Unitful.Temperature; output=[u"Pa", u"m^3", u"m^3"])
    T = float(ustrip(u"K", T))
    (P_sat, v_l, v_v) = get_sat_pure(model, T)
    return (P_sat*u"Pa", v_l*u"m^3", v_v*u"m^3")
end

function get_enthalpy_vap(model::SAFT, T::Unitful.Temperature)
    T = float(ustrip(u"K", T))
    return uconvert(output, get_enthalpy_vap(model, T)*u"J/mol")
end

function get_pressure(model::SAFT, v::Unitful.Volume, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"Pa")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    v = float(ustrip(u"m^3", v)); T = float(ustrip(u"K", T))
    return uconvert(output, get_pressure(model, z, v, T; phase=phase)*u"Pa")
end

function get_entropy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"J")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_entropy(model, p, T, z; phase=phase)*u"J")
end

function get_chemical_potential(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"J")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_chemical_potential(model, p, T, z; phase=phase)*u"J")
end

function get_internal_energy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"J")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_internal_energy(model, p, T, z; phase=phase)*u"J")
end

function get_enthalpy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"J")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_enthalpy(model, p, T, z; phase=phase)*u"J")
end

function get_Gibbs_free_energy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"J")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_Gibbs_free_energy(model, p, T, z; phase=phase)*u"J")
end

function get_Helmholtz_free_energy(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"J")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_Helmholtz_free_energy(model, p, T, z; phase=phase)*u"J")
end

function get_isochoric_heat_capacity(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"J/K")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_isochoric_heat_capacity(model, p, T, z; phase=phase)*u"J/K")
end

function get_isobaric_heat_capacity(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"J/K")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_isobaric_heat_capacity(model, p, T, z; phase=phase)*u"J/K")
end

function get_isothermal_compressibility(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"1/Pa")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? ustrip(u"mol", z) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_isothermal_compressibility(model, p, T, z; phase=phase)*u"1/Pa")
end

function get_isentropic_compressibility(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"1/Pa")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? ustrip(u"mol", z) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_isentropic_compressibility(model, p, T, z; phase=phase)*u"1/Pa")
end

function get_speed_of_sound(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"m/s")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? ustrip(u"mol", z) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_speed_of_sound(model, p, T, z; phase=phase)*u"m/s")
end

function get_isobaric_expansivity(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"1/K")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? ustrip(u"mol", z) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_isobaric_expansivity(model, p, T, z; phase=phase)*u"1/K")
end

function get_Joule_Thomson_coefficient(model::SAFT, p::Unitful.Pressure, T::Unitful.Temperature, z=[1.]; phase="unknown", output=u"K/Pa")
    z = typeof(z) <: Array{U} where {U<:Unitful.Amount} ? float(ustrip(u"mol", z)) : z
    p = float(ustrip(u"Pa", p)); T = float(ustrip(u"K", T))
    return uconvert(output, get_Joule_Thomson_coefficient(model, p, T, z; phase=phase)*u"K/Pa")
end
