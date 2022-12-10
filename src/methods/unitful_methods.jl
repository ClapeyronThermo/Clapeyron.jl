function pressure(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"Pa")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = pressure(model, _v, _T,_z)*u"Pa"
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
        function $fn(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=$unit)
            st = standarize(model,v,T,z)
            _v,_T,_z = state_to_vt(model,st)
            res = $VT_fn(model, _v, _T,_z)*$unit
            return uconvert(output, res)
        end

        function $fn(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown, output=$unit)
            st = standarize(model,p,T,z)
            _p,_T,_z = state_to_pt(model,st)
            res = $fn(model, _p, _T, _z; phase=phase)*($unit)
            return uconvert(output, res)
        end
    end
end

function volume(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown, output=u"m^3")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = volume(model, _p, _T, _z; phase=phase)*u"m^3"
    return uconvert(output, res)
end

#compressibility_factor
function compressibility_factor(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.])
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = VT_compressibility_factor(model, _v, _T,_z)
    return res
end

function compressibility_factor(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase=:unknown)
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = compressibility_factor(model, _p, _T, _z; phase=phase)
    return res
end

#second_virial_coefficient
function second_virial_coefficient(model::EoSModel, T::Unitful.Temperature, z=SA[1.];output=u"m^3")
    st = standarize(model,-1,T,z)
    _,_T,_z = state_to_pt(model,st)
    res = second_virial_coefficient(model, _T,_z)*u"m^3"
    return uconvert(output, res)
end

#pip
function pip(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.])
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = pip(model, _v, _T,_z)
    return res
end

#inversion_temperature
function inversion_temperature(model::EoSModel, p::Unitful.Pressure, z=SA[1.]; output=u"K")
    st = standarize(model,p,-1,z)
    _p,_,_z = state_to_pt(model,st)
    res = inversion_temperature(model, _p, _z)*u"K"
    return uconvert(output, res)
end

#enthalpy_vap
function enthalpy_vap(model::EoSModel, T::Unitful.Temperature; output=u"J")
    st = standarize(model,-1,T,SA[1.0])
    _,_T,_ = state_to_pt(model,st)
    res = enthalpy_vap(model, _T)*u"J"
    return uconvert(output,res)
end

function saturation_pressure(model::EoSModel, T::Unitful.Temperature; output=[u"Pa", u"m^3", u"m^3"])
    st = standarize(model,-1,T,SA[1.0])
    _,_T,_ = state_to_pt(model,st)
    (P_sat, v_l, v_v) = saturation_pressure(model,_T)
    _P_sat = uconvert(output[1],P_sat*u"Pa")
    _v_l = uconvert(output[2],v_l*u"m^3")
    _v_v = uconvert(output[3],v_v*u"m^3")
    return (_P_sat,_v_l,_v_v)
end
