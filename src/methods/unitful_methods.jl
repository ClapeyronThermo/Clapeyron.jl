function pressure(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"Pa")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = pressure(model, _v, _T,_z)*u"Pa"
    return uconvert(output, res)
end


function volume(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"m^3")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = volume(model, _p, _T, _z; phase=phase)*u"m^3"
    return uconvert(output, res)
end

#enthalpy
function enthalpy(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"J")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_enthalpy(model, _v, _T,_z)*u"J"
    return uconvert(output, res)
end


function enthalpy(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = enthalpy(model, _p, _T, _z; phase=phase)*u"J"
    return uconvert(output, res)
end

#entropy
function entropy(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"J")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_entropy(model, _v, _T,_z)*u"J"
    return uconvert(output, res)
end


function entropy(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = entropy(model, _p, _T, _z; phase=phase)*u"J"
    return uconvert(output, res)
end

#internal_energy
function internal_energy(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"J")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_internal_energy(model, _v, _T,_z)*u"J"
    return uconvert(output, res)
end


function internal_energy(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = internal_energy(model, _p, _T, _z; phase=phase)*u"J"
    return uconvert(output, res)
end

#gibbs_free_energy
function gibbs_free_energy(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"J")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_gibbs_free_energy(model, _v, _T,_z)*u"J"
    return uconvert(output, res)
end


function gibbs_free_energy(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = gibbs_free_energy(model, _p, _T, _z; phase=phase)*u"J"
    return uconvert(output, res)
end

#helmholtz_free_energy
function helmholtz_free_energy(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"J")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_helmholtz_free_energy(model, _v, _T,_z)*u"J"
    return uconvert(output, res)
end


function helmholtz_free_energy(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = helmholtz_free_energy(model, _p, _T, _z; phase=phase)*u"J"
    return uconvert(output, res)
end

#isochoric_heat_capacity
function isochoric_heat_capacity(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"J/K")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_isochoric_heat_capacity(model, _v, _T,_z)*u"J/K"
    return uconvert(output, res)
end

function isochoric_heat_capacity(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J/K")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = isochoric_heat_capacity(model, _p, _T, _z; phase=phase)*u"J/K"
    return uconvert(output, res)
end



#isobaric_heat_capacity
function isobaric_heat_capacity(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"J/K")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_isobaric_heat_capacity(model, _v, _T,_z)*u"J/K"
    return uconvert(output, res)
end

function isobaric_heat_capacity(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"J/K")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = isobaric_heat_capacity(model, _p, _T, _z; phase=phase)*u"J/K"
    return uconvert(output, res)
end

#isothermal_compressibility
function isothermal_compressibility(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"1/Pa")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_isothermal_compressibility(model, _v, _T,_z)*u"1/Pa"
    return uconvert(output, res)
end

function isothermal_compressibility(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"1/Pa")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = isothermal_compressibility(model, _p, _T, _z; phase=phase)*u"1/Pa"
    return uconvert(output, res)
end

#isentropic_compressibility
function isentropic_compressibility(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"1/Pa")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_isentropic_compressibility(model, _v, _T,_z)*u"1/Pa"
    return uconvert(output, res)
end

function isentropic_compressibility(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"1/Pa")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = isentropic_compressibility(model, _p, _T, _z; phase=phase)*u"1/Pa"
    return uconvert(output, res)
end


#speed_of_sound
function speed_of_sound(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"m/s")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_speed_of_sound(model, _v, _T,_z)*u"m/s"
    return uconvert(output, res)
end

function speed_of_sound(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"m/s")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = speed_of_sound(model, _p, _T, _z; phase=phase)*u"m/s"
    return uconvert(output, res)
end

#isobaric_expansivity
function isobaric_expansivity(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"1/K")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_isobaric_expansivity(model, _v, _T,_z)*u"1/K"
    return uconvert(output, res)
end

function isobaric_expansivity(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"1/K")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = isobaric_expansivity(model, _p, _T, _z; phase=phase)*u"1/K"
    return uconvert(output, res)
end

#joule_thomson_coefficient
function joule_thomson_coefficient(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.]; output=u"K/Pa")
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_joule_thomson_coefficient(model, _v, _T,_z)*u"K/Pa"
    return uconvert(output, res)
end

function joule_thomson_coefficient(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown", output=u"K/Pa")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = joule_thomson_coefficient(model, _p, _T, _z; phase=phase)*u"K/Pa"
    return uconvert(output, res)
end

#compressibility_factor
function compressibility_factor(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.])
    st = standarize(model,v,T,z)
    _v,_T,_z = state_to_vt(model,st)
    res = vt_joule_thomson_coefficient(model, _v, _T,_z)
    return res
end

function compressibility_factor(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    res = joule_thomson_coefficient(model, _p, _T, _z; phase=phase)
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

function sat_pure(model::EoSModel, T::Unitful.Temperature; output=[u"Pa", u"m^3", u"m^3"])
    st = standarize(model,-1,T,SA[1.0])
    _,_T,_ = state_to_pt(model,st)
    (P_sat, v_l, v_v) = sat_pure(model,_T)
    _P_sat = uconvert(output[1],P_sat*u"Pa")
    _v_l = uconvert(output[2],v_l*u"m^3")
    _v_v = uconvert(output[3],v_v*u"m^3")
    return (_P_sat,_v_l,_v_v)
end


#molar density 
function molar_density(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown",threaded=true, output=u"mol/m^3")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    
    res = molar_density(model,_p,_T,_z;phase=phase,threaded=threaded) *u"mol/m^3"
    return uconvert(output, res)
end


#mass density 
function mass_density(model::EoSModel, p::Unitful.Pressure, T::Unitful.Temperature, z=SA[1.]; phase="unknown",threaded=true, output=u"mol/m^3")
    st = standarize(model,p,T,z)
    _p,_T,_z = state_to_pt(model,st)
    
    res = mass_density(model,_p,_T,_z;phase=phase,threaded=threaded) *u"kg/m^3"
    return uconvert(output, res)
end


