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

function standardize(model,x,T,z)
    xs = standardize(x,nothing)
    ts = standardize(T,1u"K")
    zs = standardize(model,z,1u"mol")
    return StdState(xs,ts,zs)
end

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
standardize(x::Number,st::Nothing) = x #default
standardize(x::Unitful.Pressure,st::Nothing) = ustrip(u"Pa",x)
standardize(x::__VolumeKind,st::Nothing) = upreferred(x)
standardize(x::Number,st::Unitful.Temperature) = x #default
standardize(x::Unitful.Temperature,st::Unitful.Temperature) = ustrip(u"K",x)

#handle vector of compounds
standardize(model,x::Number,st::Unitful.Amount) = C.SA[x]
standardize(model,x::Unitful.Amount,st::Unitful.Amount) = C.SA[ustrip(u"mol",x)]
standardize(model,x::AbstractVector{<:Number},st::Unitful.Amount) = x
standardize(model,x::AbstractVector{<:Unitful.Amount},st::Unitful.Amount) = ustrip.(u"mol",x)
#mass forms
standardize(model,x::Unitful.Mass,st::Unitful.Amount) = C.SA[1000*ustrip(u"kg",x)/mw(model)[1]]
function standardize(model,x::AbstractVector{<:Unitful.Mass},st::Unitful.Amount)
    return map(y -> 1000*ustrip(u"kg",y[1])/y[2],zip(x,mw(model)))
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


function state_to_vt(model,st::StdState)
    V = total_volume(model,st.x,st.z)
    T = st.T
    z = st.z
    return V,T,z
end

function state_to_pt(model,st::StdState)
    p = st.x
    T = st.T
    z = st.z
    return p,T,z
end

function C.pressure(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"Pa")
    st = standardize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = C.pressure(model, _v, _T,_z)*u"Pa"
    return uconvert(output, res)
end

for (fn,unit) in [
    (:chemical_potential, u"J/mol"),
    (:chemical_potential_res, u"J/mol"),
    (:compressibility_factor, NoUnits),
    (:enthalpy, u"J"),
    (:entropy, u"J/K"),
    (:entropy_res, u"J/K"),
    (:gibbs_free_energy, u"J"),
    (:helmholtz_free_energy, u"J"),
    (:internal_energy, u"J"),
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
    @eval begin
        function C.$fn(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=$unit)
            st = standardize(model,v,T,z)
            _v,_T,_z = state_to_vt(model,st)
            res = C.$VT_fn(model, _v, _T,_z)*$unit
            return uconvert.(output, res)
        end

        function C.$fn(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown, output=$unit, threaded=true)
            st = standardize(model,p,T,z)
            _p,_T,_z = state_to_pt(model,st)
            res = C.$fn(model, _p, _T, _z; phase, threaded)*($unit)
            return uconvert.(output, res)
        end
    end
end

function C.volume(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown, output=u"m^3", threaded=true, vol0=nothing)
    st = standardize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = volume(model, _p, _T, _z; phase, threaded, vol0)*u"m^3"
    return uconvert(output, res)
end

#second_virial_coefficient
function C.second_virial_coefficient(model::EoSModel, T::Unitful.Temperature, z=SA[1.]; output=u"m^3")
    st = standardize(model,-1,T,z)
    _,_T,_z = state_to_pt(model,st)
    res = second_virial_coefficient(model, _T,_z)*u"m^3"
    return uconvert(output, res)
end

#pip
function C.pip(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.])
    st = standardize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = C.pip(model, _v, _T,_z)
    return res
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

function C.fugacity_coefficient(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown, threaded=true)
    st = standardize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = C.fugacity_coefficient(model, _p, _T, _z; phase, threaded)
    return res
end

function C.volume_virial(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; output=u"m^3")
    st = standardize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = C.volume_virial(model, _p, _T, _z)*u"m^3"
    return uconvert(output, res)
end

# x0_psat fallback method
function C.x0_psat(model::EoSModel, T::Unitful.Temperature, Tc::Unitful.Temperature, Vc::__VolumeKind; output=u"Pa")
    st = standardize(model, Vc, T, SA[1.])
    _Vc, _T, _ = state_to_vt(model,st)
    _Tc = standardize(Tc, 1u"K")
    res = C.x0_psat(model, _T, _Tc, _Vc)*u"Pa"
    return uconvert(output, res)
end

# x0_psat interface
for modeltype in (:EoSModel, :SingleFluid, :MultiFluid, :CompositeModel)
    @eval function C.x0_psat(model::($modeltype), T::Unitful.Temperature, crit=nothing; output=u"Pa")
        uconvert(output, C.x0_psat(model, standardize(T, 1u"K"), crit)*u"Pa")
    end
end

# x0_sat_pure using z
for modeltype in (:EoSModel, :SingleFluid, :ExtendedCorrespondingStates)
    @eval function C.x0_sat_pure(model::($modeltype), T::Unitful.Temperature, z=SA[1.0]; output=(u"m^3", u"m^3"))
        st = standardize(model,-1,T,z)
        _,_T,_z = state_to_pt(model,st)
        v_l, v_v = C.x0_sat_pure(model, _T, _z)
        _v_l = uconvert(output[1],v_l*u"m^3")
        _v_v = uconvert(output[2],v_v*u"m^3")
        return (_v_l,_v_v)
    end
end
# x0_sat_pure without z
for modeltype in (:CompositeModel, :MultiFluid, :LJRef, :(Clapeyron.ActivityModel), :(Clapeyron.ABCubicModel))
    @eval function C.x0_sat_pure(model::($modeltype), T::Unitful.Temperature; output=(u"m^3", u"m^3"))
        v_l, v_v = C.x0_sat_pure(model, standardize(T, 1u"K"))
        _v_l = uconvert(output[1],v_l*u"m^3")
        _v_v = uconvert(output[2],v_v*u"m^3")
        return (_v_l,_v_v)
    end
end

# x0_saturation_temperature
for modeltype in (:EoSModel, :SingleFluid, :MultiFluid)
    @eval function C.x0_saturation_temperature(model::($modeltype), p::Unitful.Pressure; output=(u"K", u"m^3", u"m^3"))
        (T_sat, v_l, v_v) = C.x0_saturation_temperature(model, standardize(p, nothing))
        _T_sat = uconvert(output[1],T_sat*u"K")
        _v_l = uconvert(output[2],v_l*u"m^3")
        _v_v = uconvert(output[3],v_v*u"m^3")
        return (_T_sat,_v_l,_v_v)
    end
end

function C.x0_volume(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown, output=u"m^3")
    st = standardize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = C.x0_volume(model, _p, _T, _z; phase)*u"m^3"
    return uconvert(output, res)
end

# x0_volume_gas
for modeltype in (:EoSModel, :MultiFluid, :SanchezLacombe, :(Clapeyron.DAPTModel))
    @eval function C.x0_volume_gas(model::($modeltype), p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; output=u"m^3")
        st = standardize(model,p,T,z)
        _p,_T,_z = state_to_pt(model,st)
        res = C.x0_volume_gas(model, _p, _T, _z)*u"m^3"
        return uconvert(output, res)
    end
end

# x0_volume_liquid with a default value for z
for modeltype in (:SingleFluid, :ExtendedCorrespondingStates, :(Clapeyron.ActivityModel), :(Clapeyron.AnalyticalSLVModel))
    @eval function C.x0_volume_liquid(model::($modeltype), T::Unitful.Temperature, z=SA[1.]; output=u"m^3")
        st = standardize(model,-1,T,z)
        _,_T,_z = state_to_pt(model,st)
        res = C.x0_volume_liquid(model, _T, _z)*u"m^3"
        return uconvert(output, res)
    end
end

# x0_volume_liquid without a default value for z
for modeltype in (:EoSModel, :MultiFluid, :SanchezLacombe, :(Clapeyron.SAFTVRQMieModel),
                  :(Clapeyron.softSAFTModel), :(Clapeyron.PeTSModel), :(Clapeyron.SAFTgammaMieModel),
                  :(Clapeyron.SAFTVRMieModel), :(Clapeyron.BACKSAFTModel))
    @eval function C.x0_volume_liquid(model::($modeltype), T::Unitful.Temperature, z; output=u"m^3")
        st = standardize(model,-1,T,z)
        _,_T,_z = state_to_pt(model,st)
        res = C.x0_volume_liquid(model, _T, _z)*u"m^3"
        return uconvert(output, res)
    end
end

# x0_volume_solid with a default value for z
function C.x0_volume_solid(model::Clapeyron.AnalyticalSLVModel, T::Unitful.Temperature, z=SA[1.]; output=u"m^3")
    st = standardize(model,-1,T,z)
    _,_T,_z = state_to_pt(model,st)
    res = C.x0_volume_solid(model, _T, _z)*u"m^3"
    return uconvert(output, res)
end

# x0_volume_solid without a default value for z
function C.x0_volume_solid(model::EoSModel, T::Unitful.Temperature, z; output=u"m^3")
    st = standardize(model,-1,T,z)
    _,_T,_z = state_to_pt(model,st)
    res = C.x0_volume_solid(model, _T, _z)*u"m^3"
    return uconvert(output, res)
end

end #module
