function _multiphase_gibbs(model,p,T,result)
    comps, β, volumes = result
    if model isa PTFlashWrapper && length(β) == 2 #TODO: in the future, PTFlashWrappers could be multiphase
        if model.model isa RestrictedEquilibriaModel
            return zero(eltype(β))
        end
    end

    res = zero(Base.promote_eltype(model,p,T,comps[1]))
    np = length(comps)
    for i in 1:np
        vi = volumes[i]
        res += β[i] * (eos(model,vi,T,comps[i]) + p*vi)
    end
    return res
end

"""
    TPFlashMethod <: ThermodynamicMethod

Abstract type for `tp_flash` routines. it requires defining `numphases(method)` and `tp_flash_impl(model,p,T,n,method)`.
If the method accept component-dependent inputs, it also should define `index_reduction(method,nonzero_indices)`

"""
abstract type TPFlashMethod <: ThermodynamicMethod end

"""
    tp_flash(model, p, T, n, method::TPFlashMethod = DETPFlash())

Routine to solve non-reactive multicomponent flash problem.
The default method uses Global Optimization. see [`DETPFlash`](@ref)

Inputs:
 - T, Temperature
 - p, Pressure
 - n, vector of number of moles of each species

Outputs - Tuple containing:
 - xᵢⱼ, Array of mole fractions of species j in phase i
 - nᵢⱼ, Array of mole numbers of species j in phase i, [mol]
 - G, Gibbs Free Energy of Equilibrium Mixture [J]
"""
function tp_flash(model::EoSModel, p, T, n; kwargs...)
    method = init_preferred_method(tp_flash,model,kwargs)
    return tp_flash(model, p, T, n, method)
end

function tp_flash(model::EoSModel, p, T, n,method::TPFlashMethod)
    result = tp_flash2(model, p, T, n,method)
    return tp_flash2_to_tpflash(model,p,T,n,result)
end

"""
    numphases(method::TPFlashMethod)

Return the number of phases supported by the TP flash method. by default its set to 2.
it the method allows it, you can set the number of phases by doing `method(;numphases = n)`.
"""
numphases(method::TPFlashMethod) = 2


"""
    supports_reduction(method::TPFlashMethod)::Bool

Checks if a TP Flash method supports index reduction (the ability to prune model components with compositions lower than a threshold). 
All current Clapeyron methods support index reduction, but some methods that alllow passing a cache could have problems
"""
supports_reduction(method::TPFlashMethod) = true
#default
include("tp_flash/DifferentialEvolutiontp_flash.jl")
include("tp_flash/michelsen_base.jl")
include("tp_flash/Michelsentp_flash.jl")
include("tp_flash/RachfordRicetp_flash.jl")
include("tp_flash/MCFlashJL.jl")
include("tp_flash/multiphase.jl")

function init_preferred_method(method::typeof(tp_flash),model::EoSModel,kwargs)
    if length(model) == 2 && any(x->haskey(kwargs,x),(:v0,:noncondensables,:nonvolatiles))
        return MichelsenTPFlash(;kwargs...)
    elseif any(x->haskey(kwargs,x),(:numphases,:max_steps,:population_size,:time_limit,:verbose,:logspace))
        return DETPFlash(;kwargs...)
    else
        return MultiPhaseTPFlash(;kwargs...)
    end
end

export tp_flash

#we use tp_flash2 and transform to tp_flash
function tp_flash2(model::EoSModel, p, T, n;kwargs...)
    method = init_preferred_method(tp_flash,model,kwargs)
    return tp_flash2(model,p,T,n,method)
end

function tp_flash2(model::EoSModel, p, T, n,method::TPFlashMethod)
    check_arraysize(model,n)
    if supports_reduction(method)
        model_r,idx_r = index_reduction(model,n)
        n_r = n[idx_r]
        method_r = index_reduction(method,idx_r)
    else
        model_r,idx_r = model,1:length(model)
        method_r,n_r = method,n
    end
    
    if length(model_r) == 1 || numphases(method) == 1
        vols = [volume(model_r,p,T,n_r)]
        β = [sum(z)*one(vols[1])]
        z = zeros(eltype(β),length(n))
        z .= n
        z ./= sum(z)
        return [z], β, vols,zero(vols)
    end
    
    comps,β,vols,g = tp_flash_impl(model_r,p,T,n_r,method_r)
    β ./= sum(β)
    β .*= sum(n)
    if supports_reduction(method)
        for i in 1:length(comps)
            xi_r = comps[i]
            xi = index_expansion(xi_r,idx_r)
            comps[i] = xi
        end
    end
    idx_sort = sortperm(vols)
    return comps[idx_sort],β[idx_sort],vols[idx_sort],g
end

function tp_flash2_to_tpflash(model,p,T,z,result)
    comps, β, volumes, g = result
    nc = length(z)
    np = length(comps)

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

