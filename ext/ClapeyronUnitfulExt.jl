module ClapeyronUnitfulExt
using Clapeyron
using Clapeyron: SA
using Unitful
const C = Clapeyron

import Unitful: @u_str

Unitful.@derived_dimension __MassDensity Unitful.𝐌/Unitful.𝐋^3
Unitful.@derived_dimension __MolDensity Unitful.𝐍/Unitful.𝐋^3
Unitful.@derived_dimension __MassVolume Unitful.𝐋^3/Unitful.𝐌
Unitful.@derived_dimension __MolVolume Unitful.𝐋^3/Unitful.𝐍
Unitful.@derived_dimension __MolAmount Unitful.𝐍

const __VolumeKind = Union{__MassDensity,__MolDensity,__MassVolume,__MolVolume,Unitful.Volume}

struct StdState{XS,TS,ZS}
    x::XS
    T::TS
    z::ZS
end

function standarize(model,x,T,z)
    xs = standarize(x,nothing)
    ts = standarize(T,1u"K")
    zs = standarize(model,z,1u"mol")
    return StdState(xs,ts,zs)
end

#utility func
_u0(x::T,val::V) where {T,V} = ustrip(uconvert(x,val))

function mw(model)
    if C.has_groups(model)
        n = length(model)
        mw_comp = zeros(n)
        v = model.groups.n_flattenedgroups
        mw_gc = C.mw(model)
        for i in eachindex(mw_comp)
            mw_comp[i] = dot(mw_gc[i],v[i])
        end
        return mw_comp
    else
        return C.mw(model)
    end
end
#handle pressure/volume/temp
standarize(x::Number,st::Nothing) = x #default
standarize(x::Unitful.Pressure,st::Nothing) = uconvert(u"Pa",x)
standarize(x::__VolumeKind,st::Nothing) = upreferred(x)
standarize(x::Number,st::Unitful.Temperature) = x #default
standarize(x::Unitful.Temperature,st::Unitful.Temperature) = _u0(u"K",x)

#handle vector of compounds
standarize(model,x::Number,st::Unitful.Amount) = C.SA[x]
standarize(model,x::Unitful.Amount,st::Unitful.Amount) = C.SA[_u0(u"mol",x)]
standarize(model,x::AbstractVector{<:Number},st::Unitful.Amount) = x
standarize(model,x::AbstractVector{<:Unitful.Amount},st::Unitful.Amount) = _u0.(u"mol",x)
#mass forms
standarize(model,x::Unitful.Mass,st::Unitful.Amount) = C.SA[1000*_u0(u"kg",x)/mw(model)[1]]
function standarize(model,x::AbstractVector{<:Unitful.Mass},st::Unitful.Amount)
    return map(y -> 1000*_u0(u"kg",y[1])/y[2],zip(x,mw(model)))
end

function mass(model,z)
    Mw = mw(model)
    res = zero(first(z)+first(Mw))
    for i in eachindex(z)
        res += 0.001*Mw[i]*z[i]
    end
    return res
end

#functions to turn any volumekind into total volume
total_volume(model,x::Number,z) = x
total_volume(model,x::Unitful.Volume,z) = ustrip(x)
total_volume(model,x::__MolVolume,z) = ustrip(x) * C.molecular_weight(model,z)
total_volume(model,x::__MolDensity,z) = C.molecular_weight(model,z) / ustrip(x)
total_volume(model,x::__MassVolume,z) = ustrip(x) * mass(model,z)
total_volume(model,x::__MassDensity,z) = mass(model,z) / ustrip(x)

_pressure(x::Number) = x
_pressure(x::Unitful.Pressure) = ustrip(x)

function state_to_vt(model,st::StdState)
    V = total_volume(model,st.x,st.z)
    T = st.T
    z = st.z
    return V,T,z
end

function state_to_pt(model,st::StdState)
    p = _pressure(st.x)
    T = st.T
    z = st.z
    return p,T,z
end

function C.pressure(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"Pa")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = C.pressure(model, _v, _T,_z)*u"Pa"
    return uconvert(output, res)
end

for (fn,unit) in Iterators.zip(
    [:enthalpy,
    :entropy,
    :internal_energy,
    :gibbs_free_energy,
    :helmholtz_free_energy,
    :isochoric_heat_capacity,
    :isobaric_heat_capacity,
    :isothermal_compressibility,
    :isentropic_compressibility,
    :speed_of_sound,
    :isobaric_expansivity,
    :joule_thomson_coefficient,
    :mass_density,
    :molar_density],
    [u"J",
    u"J/K",
    u"J",
    u"J",
    u"J",
    u"J/K",
    u"J/K",
    u"Pa^-1",
    u"Pa^-1",
    u"m/s",
    u"K^-1",
    u"K/Pa",
    u"kg/m^3",
    u"mol/m^3"])
    VT_fn = Symbol(:VT_,fn)
    @eval begin
        function C.$fn(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=$unit)
            st = standarize(model,v,T,z)
            _v,_T,_z = state_to_vt(model,st)
            res = C.$VT_fn(model, _v, _T,_z)*$unit
            return uconvert(output, res)
        end

        function C.$fn(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown, output=$unit, threaded=true)
            st = standarize(model,p,T,z)
            _p,_T,_z = state_to_pt(model,st)
            res = C.$fn(model, _p, _T, _z; phase, threaded)*($unit)
            return uconvert(output, res)
        end
    end
end

function C.volume(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown, output=u"m^3", threaded=true)
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = volume(model, _p, _T, _z; phase, threaded)*u"m^3"
    return uconvert(output, res)
end

#compressibility_factor
function C.compressibility_factor(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.])
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = C.VT_compressibility_factor(model, _v, _T,_z)
    return res
end

function C.compressibility_factor(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown, threaded=true)
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = compressibility_factor(model, _p, _T, _z; phase, threaded)
    return res
end

#second_virial_coefficient
function C.second_virial_coefficient(model::EoSModel, T::Unitful.Temperature, z=SA[1.]; output=u"m^3")
    st = standarize(model,-1,T,z)
    _,_T,_z = state_to_pt(model,st)
    res = second_virial_coefficient(model, _T,_z)*u"m^3"
    return uconvert(output, res)
end

#pip
function C.pip(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.])
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = C.pip(model, _v, _T,_z)
    return res
end

#inversion_temperature
function C.inversion_temperature(model::EoSModel, p::Unitful.Pressure, z=SA[1.]; output=u"K")
    st = standarize(model,p,-1,z)
    _p,_,_z = state_to_pt(model,st)
    res = inversion_temperature(model, _p, _z)*u"K"
    return uconvert(output, res)
end

#enthalpy_vap
function C.enthalpy_vap(model::EoSModel, T::Unitful.Temperature; output=u"J")
    st = standarize(model,-1,T,SA[1.0])
    _,_T,_ = state_to_pt(model,st)
    res = enthalpy_vap(model, _T)*u"J"
    return uconvert(output,res)
end

function C.saturation_pressure(model::EoSModel, T::Unitful.Temperature; output=(u"Pa", u"m^3", u"m^3"))
    st = standarize(model,-1,T,C.SA[1.0])
    _,_T,_ = state_to_pt(model,st)
    (P_sat, v_l, v_v) = saturation_pressure(model,_T)
    _P_sat = uconvert(output[1],P_sat*u"Pa")
    _v_l = uconvert(output[2],v_l*u"m^3")
    _v_v = uconvert(output[3],v_v*u"m^3")
    return (_P_sat,_v_l,_v_v)
end

end #module
