
function flash(m::M,model,prop1,prop2,z,args...;kwargs...) where M
    if m == PT
        result = tp_flash2(model,prop1,prop2,z,args...;kwargs...)
    elseif m == PH
        return ph_flash(model,prop1,prop2,z,args...;kwargs...)
    elseif m == PS
        return ps_flash(model,prop1,prop2,z,args...;kwargs...)
    else
        throw(error("$m flash not implemented."))
    end
end

function ph_flash(model,p,h,z,T0 = Tproperty(model,p,h,z),K0 = nothing)
    if length(model) == 1
        return ph_flash_pure(model,p,h,z,T0)
    end
    throw(error("ph_flash not implemented for multicomponent models."))
end

function ps_flash(model,p,s,z,T0 = Tproperty(model,p,s,z,entropy),K0 = nothing)
    if length(model) == 1
        return ps_flash_pure(model,p,h,z,T0)
    end
    throw(error("ps_flash not implemented for multicomponent models."))
end

function ph_flash_pure(model,p,h,z,T0 = nothing)
    T,vl,vv = saturation_temperature(model,p)
    ∑z = sum(z)
    hl = ∑z*VT_enthalpy(model,vl,T,SA[1.0])
    hv = ∑z*VT_enthalpy(model,vv,T,SA[1.0])
    βv = (h - hl)/(hv - hl)

    if !(0 <= βv <= 1)
        T,_phase = _Tproperty(model,p,h/∑z,SA[1.0],enthalpy,T0 = T0)
        v = volume(model,p,T,phase = _phase)
        return [[1.0]],[sum(z)],[v],PTFlashData(p,T,zero(v))
    end
    return build_flash_result_pure(model,p,T,z,vl,vv,βv)
end

function ps_flash_pure(model,p,s,z,T0 = nothing)
    T,vl,vv = saturation_temperature(model,p)
    ∑z = sum(z)
    sl = ∑z*VT_entropy(model,vl,T,SA[1.0])
    sv = ∑z*VT_entropy(model,vv,T,SA[1.0])
    βv = (s - sl)/(sv - sl)

    if !(0 <= βv <= 1)
        T,_phase = _Tproperty(model,p,s/∑z,SA[1.0],entropy,T0 = T0)
        v = volume(model,p,T,phase = _phase)
        return [[1.0]],[sum(z)],[v],PTFlashData(p,T,zero(v))
    end
    return build_flash_result_pure(model,p,T,z,vl,vv,βv)
end

function build_flash_result_pure(model,p,T,z,vl,vv,βv)
    comps = [[1.0],[1.0]]
    n = sum(z)
    β = [n*(1-βv),n*βv]
    volumes = [vl,vv]
    #in a single component, the gibbs energy of the (most stable) bulk phase is equal to the gibbs
    #phases of each component (only one, and equal to the bulk)
    return comps,β,volumes,PTFlashData(p,T,zero(vl))
end