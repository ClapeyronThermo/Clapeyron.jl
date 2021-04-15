function pressure(model::EoSModel, v::__VolumeKind, T::Unitful.Temperature, z=SA[1.], output=u"Pa")
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

