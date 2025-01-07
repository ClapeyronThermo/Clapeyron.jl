
function qp_f0_T!(K,z,dpdT,T0,p,T,β0)
    for i in 1:length(K)
        dlnpdTinv,logp0,T0inv = dpdT[i]
        #dTinvdlnp = -p/(dpdT[i]*T*T)
        ΔTinv = 1/T - T0inv
        K[i] = exp(ΔTinv*dlnpdTinv)
    end
    return rachfordrice(K,z) - β0
end

function qp_flash_x0(model,β,p,z,method::FlashMethod)
    if method.T0 == nothing
        if 0 <= β <= 0.01
            x = z ./ sum(z)
            T,vl,vv,y = __x0_bubble_temperature(model,p,x)
            y ./= sum(y)
            βv = β*sum(z)
            βl = sum(z) - βv
            return FlashResult(p,T,SA[x,y],SA[βl,βv],SA[vl,vv],sort = false)
        elseif 0.99 <= β <= 1.0
            y = z ./ sum(z)
            T,vl,vv,x = __x0_dew_temperature(model,p,y)
            x ./= sum(x)
            βv = β*sum(z)
            βl = sum(z) - βv
            return FlashResult(p,T,SA[x,y],SA[βl,βv],SA[vl,vv],sort = false)
        else

        pures = split_model(model)
        sat = extended_saturation_temperature.(pures,p)
        _crit = __crit_pure.(sat,pures)
        fix_sat_ti!(sat,pures,_crit,p)
        dpdT = __dlnPdTinvsat.(pures,sat,_crit,p)
        dew_prob = antoine_dew_problem(dpdT,p,z)
        Tmax = Roots.solve(dew_prob)
        bubble_prob = antoine_bubble_problem(dpdT,p,z)
        Tmin = Roots.solve(bubble_prob)
        #@show Td0,Tb0
        T0 = first.(sat)
        #Tmin,Tmax = extrema(T0)
        #we approximate sat(T) ≈ exp(-dpdT*T*T(1/T - 1/T0)/p)*p
        K = similar(T0)
        x = z ./ sum(z)
        ft(T) = qp_f0_T!(K,x,dpdT,T0,p,T,β)
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

"""
    result = qp_flash(model, q, p, n, method::FlashMethod = GeneralizedXYFlash())
    result = qp_flash(model, q, p, n; kwargs...)

Routine to solve non-reactive two-phase multicomponent flash problem. with vapour fraction - P specifications.
Wrapper around [Clapeyron.xy_flash](@ref), with automatic initial point calculations.
Inputs:
 - `q`, vapour fraction
 - `p`, pressure
 - `z`, vector of number of moles of each species

All keyword arguments are forwarded to [`GeneralizedXYFlash`](@ref).

 Outputs:
 - `result`, a [`FlashResult`](@ref) struct containing molar fractions, vapour fractions, molar volumes and the equilibrium temperature and pressure.

!!! note
    Using `qp_flash` with q = 0 or q = 1 is equivalent to calculating bubble or dew temperatures.
    Passing `GeneralizedXYFlash` as a method to [`bubble_temperature`](@ref) of [`dew_temperature`](@ref) will use `qp_flash` to calculate the bubble/dew point.
"""
function qp_flash(model::EoSModel,β,p,z;kwargs...)
    if !(0 <= β <= 1)
        throw(DomainError(β,"vapour fractions should be between 0 and 1"))
    end
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
        result1 = qflash_pure(model_r,pressure,p,β,z)
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

function bubble_temperature_impl(model::EoSModel,p,z,method::GeneralizedXYFlash)
    result = Clapeyron.qp_flash(model,0,p,z,method)
    x1,x2 = result.compositions
    v1,v2 = result.volumes
    if x1 ≈ z
        y = x2
        vl,vv = v1,v2
    else
        y = x1
        vl,vv = v2,v1
    end
    return temperature(result),vl,vv,y
end

function dew_temperature_impl(model::EoSModel,p,z,method::GeneralizedXYFlash)
    result = Clapeyron.qp_flash(model,1,p,z,method)
    x1,x2 = result.compositions
    v1,v2 = result.volumes
    if x1 ≈ z
        x = x2
        vl,vv = v2,v1
    else
        x = x1
        vl,vv = v1,v2
    end
    return temperature(result),vl,vv,x
end

export qp_flash
