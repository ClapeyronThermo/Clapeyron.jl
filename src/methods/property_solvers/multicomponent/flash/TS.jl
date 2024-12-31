"""
    result = ts_flash(model, T, s, n, method::FlashMethod = GeneralizedXYFlash())
    result = ts_flash(model, T, s, n; kwargs...)

Routine to solve non-reactive two-phase multicomponent flash problem. with T-S specifications.
Wrapper around [Clapeyron.xy_flash](@ref), with automatic initial point calculations. 
Inputs:
 - `T`, temperature
 - `s`, entropy
 - `z`, vector of number of moles of each species

All keyword arguments are forwarded to [`GeneralizedXYFlash`](@ref).

 Outputs:
 - `result`, a [`FlashResult`](@ref) struct containing molar fractions, vapour fractions, molar volumes and the equilibrium temperature and pressure.
"""
function ts_flash(model::EoSModel,T,S,z;kwargs...)
    method = init_preferred_method(ts_flash,model,kwargs)
    return ts_flash(model,T,S,z,method)
end

function init_preferred_method(method::typeof(ts_flash),model::EoSModel,kwargs) 
    GeneralizedXYFlash(;kwargs...)
end

function ts_flash(model,T,S,z,method::FlashMethod)
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
        P0 = hasfield(typeof(method),:p0) ? method.p0 : nothing
        result1 = tx_flash_pure(model,T,S,z,entropy,P0)
        return index_expansion(result1,idx_r)
    end

    result = ts_flash_impl(model_r,T,S,z_r,method_r)
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= sum(z)
    return index_expansion(result,idx_r)
end

function ts_flash_impl(model,T,S,z,method::GeneralizedXYFlash)
    flash0 = tx_flash_x0(model,T,S,z,entropy,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(entropy,S,temperature,T)
    return xy_flash(model,spec,z,flash0,method)
end

export ts_flash
