function ps_flash(model::EoSModel,p,s,z;kwargs...)
    method = init_preferred_method(ps_flash,model,kwargs)
    return ps_flash(model,p,s,z,method)
end

function init_preferred_method(method::typeof(ps_flash),model::EoSModel,kwargs) 
    GeneralizedXYFlash(;kwargs...)
end

function ps_flash(model,p,s,z,method::FlashMethod)
    check_arraysize(model,n)
    if supports_reduction(method)
        model_r,idx_r = index_reduction(model,z)
        z_r = z[idx_r]
        method_r = index_reduction(method,idx_r)
    else
        model_r,idx_r = model,trues(length(model))
        method_r,z_r = method,z
    end
    if length(model_r) == 1
        T0 = hasfield(typeof(method),:T0) ? method.T0 : nothing
        result1 = px_flash_pure(model,p,s,z,entropy,T0)
        return index_expansion(result1,idx_r)
    end
    
    result = ps_flash_impl(model_r,p,s,z_r,method_r)
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= sum(z)
    return index_expansion(result,idx_r)
end

function ps_flash_impl(model,p,s,z,method::GeneralizedXYFlash)
    flash0 = px_flash_x0(model,p,s,z,entropy,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(pressure,p,entropy,s)
    return xy_flash(model,spec,z,flash0,method)
end

export ps_flash
