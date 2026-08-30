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
C.unit_system(::Unitful.Quantity) = UnitfulJL()

function C.solve_unit(::UnitfulJL,output::Unitful.Units,output_from_f::Unitful.Units)
    uconvert(output,1*output_from_f)
end

function C.with_output_unit(res,out_units::Tuple{UnitfulJL,Unitful.Quantity})
    _,output_unit = out_units
    return res*output_unit
end

function C.with_output_unit(res,out_units::Tuple{UnitfulJL,Unitful.Units})
    _,output_unit = out_units
    return res*oneunit(output_unit)
end

C.unitful_is_pressure(::__VolumeKind) = false
C.unitful_is_pressure(::Unitful.Pressure) = true

C.uzstrip(::UnitfulJL,model,z) = __uzstrip(model,z)

__uzstrip(model,z::AbstractVector{T}) where T <: Unitful.Amount = ustrip.(u"mol",x)
__uzstrip(model,z::AbstractVector{T}) where T <: Number = x
__uzstrip(model,z::AbstractVector{T}) where T <: Unitful.Mass = map(y -> 1000*ustrip(u"kg",y[1])/y[2],zip(x,C.mw(model)))

C.uvstrip(::UnitfulJL,model,v,z) = __uvstrip(model,v,z)
__uvstrip(model,x::Number,z) = x
__uvstrip(model,x::Unitful.Volume,z) = ustrip(upreferred(x))
__uvstrip(model,x::__MolVolume,z) = ustrip(upreferred(x)) * C.molecular_weight(model,z)
__uvstrip(model,x::__MolDensity,z) = C.molecular_weight(model,z) / ustrip(x)
__uvstrip(model,x::__MassVolume,z) = ustrip(upreferred(x)) * C.molecular_weight(model,z)
__uvstrip(model,x::__MassDensity,z) = C.molecular_weight(model,z) / ustrip(upreferred(x))

#basic unit_type defs
C.unit_type(::UnitfulJL,::typeof(C.enthalpy)) = u"J"
C.unit_type(::UnitfulJL,::typeof(C.temperature)) = u"K"
C.unit_type(::UnitfulJL,::typeof(C.pressure)) = u"Pa"
C.unit_type(::UnitfulJL,::typeof(C.volume)) = u"m^3"
C.unit_type(::UnitfulJL,::typeof(C.mass_enthalpy)) = u"J"

for (fns,unit) in [
    ([:enthalpy,:enthalpy_res,:internal_energy,:internal_energy_res,:gibbs_energy,:gibbs_energy_res,:helmholtz_energy,:helmholtz_energy_res], u"J"),
    ([:mass_enthalpy,:mass_internal_energy,:mass_gibbs_energy,:mass_helmholtz_energy], u"J/kg"),
    ([:entropy,:entropy_res,:isobaric_heat_capacity,:isochoric_heat_capacity], u"J/K"),
    ([:mass_entropy,:mass_isobaric_heat_capacity,:mass_isochoric_heat_capacity], u"J/K/kg"),
    ([:isentropic_compressibility,:isothermal_compressibility], u"Pa^-1"),
    ([:isobaric_expansivity], u"K^-1"),
    ([:joule_thomson_coefficient], u"K/Pa"),
    ([:mass_density], u"kg/m^3"),
    ([:molar_density], u"mol/m^3"),
    ([:speed_of_sound], u"m/s"),
    ]
    for fn in fns
        VT_fn = Symbol(:VT_,fn)
        @eval begin
            C.unit_type(::UnitfulJL,::typeof(C.$VT_fn)) = $unit
        end
    end
end

function C.ustrip(::UnitfulJL,x::Unitful.Quantity,f::F) where F
    ustrip(uconvert(C.unit_type(UnitfulJL(),f),x))
end

#=
#inversion_temperature
function C.inversion_temperature(model::EoSModel, p::Unitful.Pressure, z=SA[1.]; output=u"K")
    st = standardize(model,p,-1,z)
    _p,_,_z = state_to_pt(model,st)
    res = inversion_temperature(model, _p, _z)*u"K"
    return uconvert(output, res)
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
end =#

end #module
