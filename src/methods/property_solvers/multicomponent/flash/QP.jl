
function qp_f0_T!(K,z,dpdT,p,T,β0)
    for i in 1:length(K)
        dlnpdTinv,logp0,T0inv = dpdT[i]
        #dTinvdlnp = -p/(dpdT[i]*T*T)
        ΔTinv = 1/T - T0inv
        K[i] = exp(ΔTinv*dlnpdTinv)
    end
    return rachfordrice(K,z) - β0
end

function qp_flash_x0(model::RestrictedEquilibriaModel,β,p,z,method::FlashMethod)
    qp_flash_x0(__tpflash_cache_model(model,p,NaN,z,:vle),β,p,z,method)
end

function qp_flash_x0(model::CompositeModel,β,p,z,method::FlashMethod)
    qp_flash_x0(model.fluid,β,p,z,method)
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
            if model isa PTFlashWrapper
                pures = model.pures
            else
                pures = split_pure_model(model)
            end
            dpdT = extended_dpdT_temperature.(pures,p)
            Tmax = antoine_dew_solve(dpdT,p,z)
            Tmin = antoine_bubble_solve(dpdT,p,z)
            #@show Td0,Tb0
            #Tmin,Tmax = extrema(T0)
            #we approximate sat(T) ≈ exp(-dpdT*T*T(1/T - 1/T0)/p)*p
            K = similar(dpdT,typeof(Tmax))
            x = z ./ sum(z)
            ft(_T) = qp_f0_T!(K,x,dpdT,p,_T,β)
            #we do a search over Tmin-Tmax domain, finding the minimum value of the objective function
            Tm = β*Tmax + (1 - β)*Tmin
            Tr1 = range(Tmin,Tm,5*length(model))
            Tr2 = range(Tm,Tmax,5*length(model))
            δβ1_min,i1 = findmin(w -> abs(ft(w)),Tr1)
            δβ2_min,i2 = findmin(w -> abs(ft(w)),Tr2)
            T00 = δβ1_min < δβ2_min ? Tr1[i1] : Tr2[i2]
            prob = Roots.ZeroProblem(ft,T00)
            T = Roots.solve(prob)
        end
    else
        T = method.T0
    end
    update_temperature!(model,T)
    r = pt_flash_x0(model,p,T,z,method)
    return r
end

"""
    result = qp_flash(model, q, p, n, method::FlashMethod = GeneralizedXYFlash())
    result = qp_flash(model, q, p, n; kwargs...)

Routine to solve non-reactive two-phase multicomponent flash problem. with vapour fraction - P specifications.
Wrapper around [Clapeyron.xy_flash](@ref), with automatic initial point calculations.
Inputs:
 - `q`, vapour fraction
 - `p`, pressure `[Pa]`
 - `n`, vector of number of moles of each species `[mol]`

All keyword arguments are forwarded to [`GeneralizedXYFlash`](@ref).

 Outputs:
 - `result`, a [`FlashResult`](@ref) struct containing molar fractions, vapour fractions, molar volumes and the equilibrium temperature and pressure.

!!! note
    Using `qp_flash` with q = 0 or q = 1 is equivalent to calculating bubble or dew temperatures.
    Passing `GeneralizedXYFlash` as a method to [`bubble_temperature`](@ref) of [`dew_temperature`](@ref) will use `qp_flash` to calculate the bubble/dew point.
"""
function qp_flash(model::EoSModel,β,p,z = SA[1.0];kwargs...)
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

    if z isa SingleComp || length(model) == 1
        result1 = qflash_pure(model,pressure,p,β,z)
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
        result1r = qflash_pure(model_r,pressure,p,β,z1r)
        return index_expansion(result1r,idx_r)
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

function qp_to_bubbledew(model,p,z,method,bubble)
    β = bubble ? 0 : 1
    result = Clapeyron.qp_flash(model,β,p,z,method)
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

qp_to_dewp(model,p,z,method) = qp_to_bubbledew(model,p,z,method,false)
qp_to_bubblep(model,p,z,method) = qp_to_bubbledew(model,p,z,method,true)


function bubble_temperature_impl(model::EoSModel,p,z,method::GeneralizedXYFlash)
    return qp_to_bubblep(model,p,z,method)
end

function dew_temperature_impl(model::EoSModel,p,z,method::GeneralizedXYFlash)
    return qp_to_dewp(model,p,z,method)
end

include("../tp_flash/RRQXFlash.jl")

export qp_flash
