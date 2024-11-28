
function qp_f0_T!(K,z,dpdT,T0,p,T,β0)
    for i in 1:length(K)
        dTinvdlnp = -p/(dpdT[i]*T*T)
        ΔTinv = 1/T - 1/T0[i]
        K[i] = exp(ΔTinv/dTinvdlnp)
    end
    return rachfordrice(K,z) - β0
end

function qp_flash_x0(model,β,p,z,method::FlashMethod)
    if method.T0 == nothing
        pures = split_model(model)
        sat = extended_saturation_temperature.(pures,p)
        
        T0 = first.(sat)
        Tmin,Tmax = extrema(T0)
        if β == 1
            T = Tmax
        elseif β == 0
            T = Tmin
        else
            #we approximate sat(T) ≈ exp(-dpdT*T*T(1/T - 1/T0)/p)*p
            dpdT = map(i -> dpdT_pure(pures[i],sat[i][3],sat[i][2],sat[i][1]),1:length(model))
            K = similar(dpdT)
            ft(T) = qp_f0_T!(K,z,dpdT,T0,p,T,β)     
            #we do a search over Tmin-Tmax domain, finding the minimum value of the objective function
            Tm = β*Tmax + (1 - β)*Tmin
            Tr1 = range(Tmin,Tm,5*length(model))
            Tr2 = range(Tm,Tmax,5*length(model))
            δβ1 = abs.(ft.(Tr1))
            δβ2 = abs.(ft.(Tr2))
                δβ1_min,i1 = findmin(δβ1)
                δβ2_min,i2 = findmin(δβ2)
                if δβ1_min < δβ2_min
                    T00 = Tr1[i1]
                else
                    T00 = Tr2[i2]
                end
                prob = Roots.ZeroProblem(ft,T00)
                T = Roots.solve(prob)
            
        end
    else
        T = method.T0
    end
    r = pt_flash_x0(model,p,T,z,method;k0 = :suggest)
    return r
end

function qp_flash(model::EoSModel,β,p,z;kwargs...)
    method = init_preferred_method(qp_flash,model,kwargs)
    return qp_flash(model,β,p,z,method)
end

function init_preferred_method(method::typeof(qp_flash),model::EoSModel,kwargs)
    GeneralizedXYFlash(;kwargs...)
end

function qp_flash(model,β,p,z,method::FlashMethod)
    check_arraysize(model,z)
    if supports_reduction(method)
        model_r,idx_r = index_reduction(model,z)
        z_r = z[idx_r]
        method_r = index_reduction(method,idx_r)
    else
        model_r,idx_r = model,trues(length(model))
        method_r,z_r = method,z
    end
    if length(model_r) == 1
        result1 = βflash_pure(model_r,pressure,p,βv,z)
        return index_expansion(result1,idx_r)
    end

    result = qp_flash_impl(model_r,β,p,z_r,method_r)
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= sum(z)
    return index_expansion(result,idx_r)
end

function qp_flash_impl(model,β,p,z,method::GeneralizedXYFlash)
    flash0 = qp_flash_x0(model,β,p,z,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(Vfrac(2),β,pressure,p)
    return xy_flash(model,spec,z,flash0,method)
end
