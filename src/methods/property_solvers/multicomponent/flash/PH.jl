function ph_flash(model,p,h,z,T0 = nothing,K0 = nothing)
    if length(model) == 1
        return ph_flash_pure(model,p,h,z,T0)
    end

    if T0 == nothing
        T,_phase = _Tproperty(model,p,h,z,enthalpy)
    else
        T = T0
        _phase = :eq #we suppose this
    end

    TT = Base.promote_eltype(model,p,h,z,T)
    if _phase != :eq 
        xx = Vector{TT}(undef,length(z))
        xx .= z
        ∑z = sum(z)
        xx ./= ∑z
        return [xx],[∑z],[volume(model,p,T,z)],PTFlashData(promote(p,T,zero(∑z))...)
    end

    if K0 == nothing
        K = suggest_K(model,p,T,z)
    else
        K = K0
    end

    spec = FlashSpecifications(pressure,p,enthalpy,h)
    βv,singlephase,_ = rachfordrice_β0(K,z)
    βv = βv*one(TT)
    #if singlephase == true, maybe initial K values overshoot the actual phase split.
    if singlephase
        Kmin,Kmax = extrema(K)
        if !(Kmin >= 1 || Kmax <= 1)
            #valid K, still single phase.
            g0 = dot(z, K) - 1. #rachford rice, supposing β = 0
            g1 = 1. - sum(zi/Ki for (zi,Ki) in zip(z,K)) #rachford rice, supposing β = 1
            if g0 <= 0 && g1 < 0 #bubble point.
                βv = eps(typeof(βv))
                singlephase = false
            elseif g0 > 0 && g1 >= 0 #dew point
                βv = one(βv) - eps(typeof(βv))
                singlephase = false
            end
        end
    else
        βv = rachfordrice(K, z; β0=βv)
    end
    y = rr_flash_vapor(K,z,βv)
    y ./= sum(y)
    x = rr_flash_liquid(K,z,βv)
    x ./= sum(x)
    βl = 1 - βv
    comps0 = [x,y]
    volumes0 = [volume(model,p,T,x,phase = :l),volume(model,p,T,y,phase = :v)]
    β0 = [βl,βv]
    return __xy_flash(model,spec,z,comps0,β0,volumes0,T)
end

function ph_flash_pure(model,p,h,z,T0 = nothing)
    Ts,vl,vv = saturation_temperature(model,p)
    ∑z = sum(z)
    hl = ∑z*VT_enthalpy(model,vl,T,SA[1.0])
    hv = ∑z*VT_enthalpy(model,vv,T,SA[1.0])
    βv = (h - hl)/(hv - hl)
    if (0 <= βv <= 1)
    elseif βv > 1
        return build_flash_result_pure(model,p,Ts,z,vl,vv,βv)
    elseif !isfinite(βv)
        [[βv]],[sum(βv)],[βv],PTFlashData(p,T,zero(v))
    elseif βv < 0 || βv > 1
        phase0 = βv < 0 ? :liquid : :vapour
        T,_phase = _Tproperty(model,p,h/∑z,SA[1.0],enthalpy,T0 = T0,phase = phase0)
        v = volume(model,p,T,phase = _phase)
        return [[1.0]],[sum(z)],[v],PTFlashData(p,T,zero(v))
    else
        build_flash_result_pure(model,p,Ts,z,vl,vv,βv)
    end
end
