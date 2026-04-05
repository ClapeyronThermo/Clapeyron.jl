"""
    result = vt_flash(model, v, T, n, method::FlashMethod = GeneralizedXYFlash())
    result = vt_flash(model, v, T, n; kwargs...)

Routine to solve non-reactive two-phase multicomponent flash problem. with V-T specifications.
Wrapper around [Clapeyron.xy_flash](@ref), with automatic initial point calculations. 
Inputs:
 - `v`, volume `[m³]`
 - `T`, temperature `[K]`
 - `n`, vector of number of moles of each species `[mol]`

All keyword arguments are forwarded to [`GeneralizedXYFlash`](@ref).

 Outputs:
 - `result`, a [`FlashResult`](@ref) struct containing molar fractions, vapour fractions, molar volumes and the equilibrium temperature and pressure.
"""
function vt_flash(model::EoSModel,V,T,z = SA[1.0];kwargs...)
    method = init_preferred_method(vt_flash,model,kwargs)
    return vt_flash(model,V,T,z,method)
end

function init_preferred_method(method::typeof(vt_flash),model::EoSModel,kwargs) 
    GeneralizedXYFlash(;kwargs...)
end

function vt_flash(model,V,T,z,method::FlashMethod)
    check_arraysize(model,z)

    if z isa SingleComp || length(model) == 1
        z1 = SVector(z[1])
        P0 = hasfield(typeof(method),:p0) ? method.p0 : nothing
        result1r = tx_flash_pure(model,T,V,z1,volume,P0)
        return result1r
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
        P0 = hasfield(typeof(method),:p0) ? method.p0 : nothing
        result1 = tx_flash_pure(model_r,T,V,z1r,volume,P0)
        return index_expansion(result1,idx_r)
    end
    
    result = vt_flash_impl(model_r,V,T,z_r,method_r)
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= sum(z)
    return index_expansion(result,idx_r)
end

function vt_flash_x0(model,V,T,z,method::GeneralizedXYFlash)
    flash0 = tx_flash_x0(model,T,V,z,volume,method)
    xl,xv = flash0.compositions
    vl,vv = flash0.volumes
    nl,nv = flash0.fractions
    n = sum(z)
    P0 = xv ./ xl .* flash0.data.p
    vbulk = V/n
    K = similar(P0)
    RT = Rgas(model)*T
    #=
    vl + vv = V
    vv = V - vl(p)
    v(p) = ng*RT/p = V - Vl(p)
    
    (vz - (1 - b)*vl)/RT = b/p
    (vz - vl + b*vl)/RT = b/p
    vz/RT - vl/RT + b*vl/RT - b/p = 0
    (vz - vl)/RT + b(vl/RT - 1/p) = 0
    (vl - vz)/(vl - RT/p) = b
    (vz - vl)(RT/p - vl) = b
    
    =#
    function f(p)
        b = (vbulk - vl)/(RT/p - vl)
        (1 - b)*vl + b*RT/p - vbulk
        K .= P0 ./ p
        return rachfordrice(K,z) - b
    end

    prob = Roots.ZeroProblem(f,flash0.data.p)

    pp = Roots.solve(prob)
    K .= P0 ./ pp

    if method.verbose
        @info "VT-flash: supposing ideal gas + constant liquid, p = $pp"
    end
    b = rachfordrice(K,z)
    rr_flash_vapor!(xv,K,z,b)
    rr_flash_liquid!(xl,K,z,b)
    fracs = SVector(n*(1-b),n*b)
    volumes = SVector(vl,(vbulk - (1-b)*vl)/b)
    return FlashResult(flash0.compositions,fracs,volumes,FlashData(pp,T))
end

function vt_flash_impl(model,V,T,z,method::GeneralizedXYFlash)
    flash0 = vt_flash_x0(model,V,T,z,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(volume,V,temperature,T)
    return xy_flash(model,spec,z,flash0,method)
end

export vt_flash
