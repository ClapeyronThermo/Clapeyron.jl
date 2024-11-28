function qt_f0_p!(K,z,p,ps,β0)
    K .= ps ./ p
    return rachfordrice(K,z) - β0
end


function qt_flash_x0(model,β,T,z,method::FlashMethod)
    if method.p0 == nothing
        pures = split_model(model)
        sat = extended_saturation_pressure.(pures,T)
        ps = first.(sat)
        K = similar(ps)
        pmin,pmax = extrema(ps)
        p0 = β*pmin + (1-β)*pmax
        fp(p) = qt_f0_p!(K,z,p,ps,β)
        prob = Roots.ZeroProblem(fp,p0)
        p = Roots.solve(prob)
    else
        p = method.p0
    end
    res =  pt_flash_x0(model,p,T,z,method;k0 = :suggest)
    return res
end

function qt_flash(model::EoSModel,β,T,z;kwargs...)
    method = init_preferred_method(qt_flash,model,kwargs)
    return qt_flash(model,β,T,z,method)
end

function init_preferred_method(method::typeof(qt_flash),model::EoSModel,kwargs) 
    GeneralizedXYFlash(;kwargs...)
end

function qt_flash(model,β,T,z,method::FlashMethod)
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
        result1 = βflash_pure(model,temperature,T,βv,z)
        return index_expansion(result1,idx_r)
    end
    
    result = qt_flash_impl(model_r,β,T,z_r,method_r)
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= sum(z)
    return index_expansion(result,idx_r)
end

function qt_flash_impl(model,β,T,z,method::GeneralizedXYFlash)
    flash0 = qt_flash_x0(model,β,T,z,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(Vfrac(2),β,temperature,T)
    return xy_flash(model,spec,z,flash0,method)
end