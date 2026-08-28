module ClapeyronUnitfulExt
using Clapeyron
using Clapeyron: SA
using Unitful
const C = Clapeyron

import Unitful: @u_str

struct UnitfulJL end

Unitful.@derived_dimension __MassDensity Unitful.𝐌/Unitful.𝐋^3
Unitful.@derived_dimension __MolDensity Unitful.𝐍/Unitful.𝐋^3
Unitful.@derived_dimension __MassVolume Unitful.𝐋^3/Unitful.𝐌
Unitful.@derived_dimension __MolVolume Unitful.𝐋^3/Unitful.𝐍
Unitful.@derived_dimension __MolAmount Unitful.𝐍

const __VolumeKind = Union{__MassDensity,__MolDensity,__MassVolume,__MolVolume,Unitful.Volume}

C.unit_system(::Unitful.Units) = UnitfulJL()
C.solve_unit(::UnitfulJL,output,output_from_f) = uconvert(output,oneunit(output_from_f))

function C.with_output_unit(res,out_units::Tuple{UnitfulJL,Unitful.Quantity})
    _,output_unit = out_units
    return res*output_unit
end

function C.with_output_unit(res,out_units::Tuple{UnitfulJL,Unitful.Units})
    _,output_unit = out_units
    return res*oneunit(output_unit)
end

C.ustrip(::UnitfulJL,x::Unitful.Quantity,typeof(C.pressure)) = ustrip(u"Pa",x)
C.ustrip(::UnitfulJL,x::Unitful.Quantity,typeof(C.temperature)) = ustrip(u"K",x)
C.uzstrip(::UnitfulJL,model,z) = __uzstrip(model,z)

__uzstrip(model,z::AbstractVector{T}) where T <: Unitful.Amount = ustrip.(u"mol",x)
__uzstrip(model,z::AbstractVector{T}) where T <: Number = x
__uzstrip(model,z::AbstractVector{T}) where T <: Unitful.Mass = map(y -> 1000*ustrip(u"kg",y[1])/y[2],zip(x,mw(model)))

C.uvstrip(::UnitfulJL,model,v,z) = __uvstrip(model,v,z)
__uvstrip(model,x::Number,z) = x
__uvstrip(model,x::Unitful.Volume,z) = ustrip(x)
__uvstrip(model,x::__MolVolume,z) = ustrip(x) * C.molecular_weight(model,z)
__uvstrip(model,x::__MolDensity,z) = C.molecular_weight(model,z) / ustrip(x)
__uvstrip(model,x::__MassVolume,z) = ustrip(x) * C.molecular_weight(model,z)
__uvstrip(model,x::__MassDensity,z) = C.molecular_weight(model,z) / ustrip(x)

for (fns,unit) in [
    ([:chemical_potential,:chemical_potential_res], u"J/mol"),
    ([:enthalpy,:enthalpy_res,:internal_energy,:internal_energy_res,:gibbs_free_energy,:gibbs_free_energy_res,:helmholtz_free_energy,:helmholtz_free_energy_res], u"J"),
    ([:enthalpy,:internal_energy,:gibbs_free_energy,:helmholtz_free_energy], u"J"),
    ([:entropy,:entropy_res], u"J/K"),

    (:isentropic_compressibility, u"Pa^-1"),
    (:isobaric_expansivity, u"K^-1"),
    (:isobaric_heat_capacity, u"J/K"),
    (:isochoric_heat_capacity, u"J/K"),
    (:isothermal_compressibility, u"Pa^-1"),
    (:joule_thomson_coefficient, u"K/Pa"),
    (:mass_density, u"kg/m^3"),
    (:molar_density, u"mol/m^3"),
    (:speed_of_sound, u"m/s"),
    ]
    VT_fn = Symbol(:VT_,fn)
end


#inversion_temperature
function C.inversion_temperature(model::EoSModel, p::Unitful.Pressure, z=SA[1.]; output=u"K")
    st = standardize(model,p,-1,z)
    _p,_,_z = state_to_pt(model,st)
    res = inversion_temperature(model, _p, _z)*u"K"
    return uconvert(output, res)
end

#enthalpy_vap
function C.enthalpy_vap(model::EoSModel, T::Unitful.Temperature; output=u"J")
    st = standardize(model,-1,T,SA[1.0])
    _,_T,_ = state_to_pt(model,st)
    res = enthalpy_vap(model, _T)*u"J"
    return uconvert(output,res)
end

function C.saturation_pressure(model::EoSModel, T::Unitful.Temperature; output=(u"Pa", u"m^3", u"m^3"))
    st = standardize(model,-1,T,C.SA[1.0])
    _,_T,_ = state_to_pt(model,st)
    (P_sat, v_l, v_v) = saturation_pressure(model,_T)
    _P_sat = uconvert(output[1],P_sat*u"Pa")
    _v_l = uconvert(output[2],v_l*u"m^3")
    _v_v = uconvert(output[3],v_v*u"m^3")
    return (_P_sat,_v_l,_v_v)
end

function C.saturation_temperature(model::EoSModel, p::Unitful.Pressure; output=(u"K", u"m^3", u"m^3"))
    (T_sat, v_l, v_v) = C.saturation_temperature(model, standardize(p, nothing))
    _T_sat = uconvert(output[1],T_sat*u"K")
    _v_l = uconvert(output[2],v_l*u"m^3")
    _v_v = uconvert(output[3],v_v*u"m^3")
    return (_T_sat,_v_l,_v_v)
end

end #module
