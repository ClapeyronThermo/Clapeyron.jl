function ps_flash(model,p,s,z,T0 = Tproperty(model,p,s,z,entropy),K0 = nothing)
    if length(model) == 1
        return ps_flash_pure(model,p,h,z,T0)
    end
    throw(error("ps_flash not implemented for multicomponent models."))
end


function ps_flash_pure(model,p,s,z,T0 = nothing)
    T,vl,vv = saturation_temperature(model,p)
    ∑z = sum(z)
    sl = ∑z*VT_entropy(model,vl,T,SA[1.0])
    sv = ∑z*VT_entropy(model,vv,T,SA[1.0])
    βv = (s - sl)/(sv - sl)

    if !(0 <= βv <= 1)
        T,_phase = _Tproperty(model,p,s/∑z,SA[1.0],entropy,T0 = T0)
        return FlashResult(model,p,T,only(z),phase = _phase)
    end
    return build_flash_result_pure(model,p,T,z,vl,vv,βv)
end
