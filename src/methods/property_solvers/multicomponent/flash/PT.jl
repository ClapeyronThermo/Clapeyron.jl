"""
    TPFlashMethod <: ThermodynamicMethod

Abstract type for `tp_flash` routines.
"""
abstract type TPFlashMethod <: FlashMethod end

"""
    tp_flash(model, p, T, n, method::TPFlashMethod = DETPFlash())

Routine to solve non-reactive multicomponent flash problem.
The default method uses Global Optimization. see [`DETPFlash`](@ref)

Inputs:
 - T, Temperature `[K]`
 - p, Pressure `[Pa]`
 - n, vector of number of moles of each species `[mol]`

Outputs - Tuple containing:
 - xᵢⱼ, Array of mole fractions of species j in phase i
 - nᵢⱼ, Array of mole numbers of species j in phase i, `[mol]`
 - G, Gibbs energy of Equilibrium Mixture `[J]`
"""
function tp_flash(model::EoSModel, p, T, n = SA[1.0]; kwargs...)
    method = init_preferred_method(tp_flash,model,kwargs)
    return tp_flash(model, p, T, n, method)
end

function tp_flash(model::EoSModel, p, T, n, method::FlashMethod)
    result = tp_flash2(model, p, T, n,method)
    return tp_flash2_to_tpflash(model,p,T,n,result)
end

#default
include("../tp_flash/DifferentialEvolutiontp_flash.jl")
include("../tp_flash/michelsen_base.jl")
include("../tp_flash/Michelsentp_flash.jl")
include("../tp_flash/RachfordRicetp_flash.jl")
include("../tp_flash/MCFlashJL.jl")
include("../tp_flash/multiphase.jl")

function init_preferred_method(method::typeof(tp_flash),model::EoSModel,kwargs) 
    if length(kwargs) == 0
        return MichelsenTPFlash()
    elseif any(x->haskey(kwargs,x),(:v0,:noncondensables,:nonvolatiles,:x0,:y0,:K0,:equilibrium))
        return MichelsenTPFlash(;kwargs...)
    elseif any(x->haskey(kwargs,x),(:numphases,:max_steps,:population_size,:time_limit,:verbose,:logspace))
        return DETPFlash(;kwargs...)
    elseif any(x->haskey(kwargs,x),(:n0,:full_tpd,:max_phases,:phase_iters))
        return MultiPhaseTPFlash(;kwargs...)
    else
        MultiPhaseTPFlash(;kwargs...)
    end
end

export tp_flash

#we use tp_flash2 and transform to tp_flash
function tp_flash2(model::EoSModel, p, T, n; kwargs...)
    method = init_preferred_method(tp_flash,model,kwargs)
    return tp_flash2(model,p,T,n,method)
end

function tp_flash2(model::EoSModel, p, T, n,method::FlashMethod)
    check_arraysize(model,n)
    
    if n isa SingleComp || length(model) == 1
        return FlashResult(model,p,T,SVector(n[1]))
    end
    
    if supports_reduction(method)
        model_r,idx_r = index_reduction(model,n)
        n_r = n[idx_r]
        method_r = index_reduction(method,idx_r)
    else
        model_r,idx_r = model,trues(length(model))
        method_r,n_r = method,n
    end
    
    if length(model_r) == 1 || numphases(method) == 1
        return FlashResult(model_r,p,T,n_r)
    end
    ∑n = sum(n_r)
    z_r = n_r ./ ∑n
    if has_a_res(model)
        λmodel,λp,λT,λz = primalval(model_r),primalval(p),primalval(T),primalval(z_r)
        λresult = tp_flash_impl(λmodel,λp,λT,λz,method_r)
        tup = (model,p,T,z)
        λtup = (λmodel,λp,λT,λz)
        result = tp_flash_ad(λresult,tup,λtup)
    else
        result = tp_flash_impl(model_r,p,T,z_r,method_r)
    end
    if !issorted(result.volumes)
        #this is in case we catch a bad result.
        result = FlashResult(result)
    end
    ∑β = sum(result.fractions)
    result.fractions ./= ∑β
    result.fractions .*= ∑n
    return index_expansion(result,idx_r)
end

function tp_flash2_to_tpflash(model,p,T,z,result)
    comps, β, volumes, data = result
    nc = length(z)
    np = length(comps)
    g = data.g
    x = similar(comps[1],(np,nc))
    n = similar(x)
    for i in 1:np
        xi = comps[i]
        βi = β[i]
        for j in 1:nc
            x[i,j] = xi[j]
            n[i,j] = xi[j]*βi
        end
    end
    return x,n,g
end

function tp_flash_impl(model,p,T,z,method::GeneralizedXYFlash)
    flash0 = px_flash_x0(model,p,T,z,temperature,method)
    isone(numphases(flash0)) && return flash0
    spec = FlashSpecifications(pressure,p,temperature,T)
    return xy_flash(model,spec,z,flash0,method)
end
