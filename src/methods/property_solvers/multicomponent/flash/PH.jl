
"""
    result = ph_flash(model, p, h, n, method::FlashMethod = GeneralizedXYFlash())
    result = ph_flash(model, p, h, n; kwargs...)

Routine to solve non-reactive two-phase multicomponent flash problem. with P-H specifications.
Wrapper around [Clapeyron.xy_flash](@ref), with automatic initial point calculations.
Inputs:
 - `p`, pressure `[Pa]`
 - `h`, enthalpy `[J]`
 - `n`, vector of number of moles of each species `[mol]`

All keyword arguments are forwarded to [`GeneralizedXYFlash`](@ref).

 Outputs:
 - `result`, a [`FlashResult`](@ref) struct containing molar fractions, vapour fractions, molar volumes and the equilibrium temperature and pressure.
"""
function ph_flash end

function ph_flash(model::EoSModel,p,h,z =SA[1.0];kwargs...)
    method = init_preferred_method(ph_flash,model,kwargs)
    return ph_flash(model,p,h,z,method)
end

function init_preferred_method(method::typeof(ph_flash),model::EoSModel,kwargs)
    GeneralizedXYFlash(;kwargs...)
end

function ph_flash(model,p,h,z,method::FlashMethod)
    check_arraysize(model,z)

    if z isa SingleComp || length(model) == 1
        z1 = SVector(z[1])
        T0 = hasfield(typeof(method),:T0) ? method.T0 : nothing
        result1 = px_flash_pure(model,p,h,z1,enthalpy,T0)
        return result1
    end

    if supports_reduction(method)
        model_r,idx_r = index_reduction(model,z)
        z_r = z[idx_r]
        method_r = index_reduction(method,idx_r)
    else
        model_r,idx_r = model,trues(length(model))
        method_r,z_r = method,z
    end

    if length(model_r) == 1
        z1r = SVector(z_r[1])
        T0 = hasfield(typeof(method),:T0) ? method.T0 : nothing
        result1r = px_flash_pure(model_r,p,h,z1r,enthalpy,T0)
        return index_expansion(result1r,idx_r)
    end

    #result = ph_flash_impl(model_r,p,h,z_r,method_r)
    if has_a_res(model)
        λmodel,λp,λh,λz = primalval(model_r),primalval(p),primalval(h),primalval(z_r)
        λresult = ph_flash_impl(λmodel,λp,λh,λz,primalval(method_r))
        tup = (model_r,p,h,z_r)
        λtup = (λmodel,λp,λh,λz)
        result = xy_flash_ad(λresult,tup,λtup,pressure,enthalpy)
    else
        result = ph_flash_impl(model_r,p,h,z,method_r)
    end
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= sum(z)
    return index_expansion(result,idx_r)
end

function ph_flash_impl(model,p,h,z,method::GeneralizedXYFlash)
    flash0 = px_flash_x0(model,p,h,z,enthalpy,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(pressure,p,enthalpy,h)
    return xy_flash(model,spec,z,flash0,method)
end

export ph_flash
