function qt_f0_p!(K,z,p,ps,β0)
    K .= ps ./ p
    return rachfordrice(K,z) - β0
end

function qt_flash_x0(model,β,T,z,method::FlashMethod)
    if method.p0 == nothing
        if 0 <= β <= 0.01
            x = z ./ sum(z)
            p,vl,vv,y = __x0_bubble_pressure(model,T,x)
            y ./= sum(y)

            return FlashResult(p,T,SA[x,y],SA[1.0-β,1.0*β],SA[vl,vv],sort = false)
        elseif 0.99 <= β <= 1.0
            y = z ./ sum(z)
            p,vl,vv,x = __x0_dew_pressure(model,T,y)
            x ./= sum(x)
            return FlashResult(p,T,SA[x,y],SA[1.0-β,1.0*β],SA[vl,vv],sort = false)
        else
            pures = split_model(model)
            sat = extended_saturation_pressure.(pures,T)
            ps = first.(sat)
            K = similar(ps)
            pmin,pmax = extrema(ps)
            fp(p) = qt_f0_p!(K,z,p,ps,β)
            pm = β*pmin + (1-β)*pmax
            pr1 = range(pmin,pm,5*length(model))
            pr2 = range(pm,pmax,5*length(model))
            δβ1 = abs.(fp.(pr1))
            δβ2 = abs.(fp.(pr2))
            δβ1_min,i1 = findmin(δβ1)
            δβ2_min,i2 = findmin(δβ2)
            if δβ1_min < δβ2_min
                p00 = pr1[i1]
            else
                p00 = pr2[i2]
            end
            prob = Roots.ZeroProblem(fp,p00)
            p = Roots.solve(prob)
        end
    else
        p = method.p0
    end
    res =  pt_flash_x0(model,p,T,z,method;k0 = :suggest)
    return res
end

"""
    result = qt_flash(model, q, T, n, method::FlashMethod = GeneralizedXYFlash())
    result = qt_flash(model, q, T, n; kwargs...)

Routine to solve non-reactive two-phase multicomponent flash problem. with vapour fraction - T specifications.
Wrapper around [Clapeyron.xy_flash](@ref), with automatic initial point calculations. 
Inputs:
 - `q`, vapour fraction
 - `T`, temperature
 - `z`, vector of number of moles of each species

All keyword arguments are forwarded to [`GeneralizedXYFlash`](@ref).

 Outputs:
 - `result`, a [`FlashResult`](@ref) struct containing molar fractions, vapour fractions, molar volumes and the equilibrium temperature and pressure.

!!! note
    Using `qt_flash` with q = 0 or q = 1 is equivalent to calculating bubble or dew pressures.
    Passing `GeneralizedXYFlash` as a method to [`bubble_pressure`](@ref) of [`dew_pressure`](@ref) will use `qt_flash` to calculate the bubble/dew point.
"""
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
        result1 = qflash_pure(model,temperature,T,βv,z)
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

function bubble_pressure_impl(model::EoSModel,T,z,method::GeneralizedXYFlash)
    result = Clapeyron.qt_flash(model,0,T,z,method)
    x1,x2 = result.compositions
    v1,v2 = result.volumes
    if x1 ≈ z
        y = x2
        vl,vv = v1,v2
    else
        y = x1
        vl,vv = v2,v1
    end
    return pressure(result),vl,vv,y
end

function dew_pressure_impl(model::EoSModel,T,z,method::GeneralizedXYFlash)
    result = Clapeyron.qt_flash(model,1,T,z,method)
    x1,x2 = result.compositions
    v1,v2 = result.volumes
    if x1 ≈ z
        x = x2
        vl,vv = v2,v1
    else
        x = x1
        vl,vv = v1,v2
    end
    return pressure(result),vl,vv,x
end

export qt_flash
