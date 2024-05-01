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
function tp_flash(model::EoSModel, p, T, n;kwargs...)
    method = init_preferred_method(tp_flash,model,kwargs)
    return tp_flash(model, p, T, n,method)
end

function tp_flash(model::EoSModel, p, T, n,method::TPFlashMethod)
    numspecies = length(model)
    if numspecies != length(n)
        error("There are ", numspecies,
            " species in the model, but the number of mole numbers specified is ",
            length(n))
    end
    if supports_reduction(method)
        model_r,idx_r = index_reduction(model,n)
        n_r = n[idx_r]
        method_r = index_reduction(method,idx_r)
    else
        model_r,idx_r = model,1:length(model)
        method_r,n_r = method,n
    end
    
    if length(model_r) == 1
        V = volume(model_r,p,T,n_r)
        return (n, n / sum(n), VT_gibbs_free_energy(model_r, V, T, n_r))
    end
    if numphases(method) == 1
        V = volume(model_r,p,T,n_r,phase =:stable)
        return (n, n / sum(n), VT_gibbs_free_energy(model_r, V, T, n_r))
    end
    
    xij_r,nij_r,g = tp_flash_impl(model_r,p,T,n_r,method_r)
    #TODO: perform stability check ritht here:
    #expand reduced model:
    if supports_reduction(method)
        nij = index_expansion(nij_r,idx_r)
        xij = index_expansion(xij_r,idx_r)
        return xij,nij,g
    else
        return xij_r,nij_r,g
    end
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
include("tp_flash/michelsen2.jl")

function init_preferred_method(method::typeof(tp_flash),model::EoSModel,kwargs)
    if haskey(kwargs,:equilibrium) || haskey(kwargs,:K0) || haskey(kwargs,:y0) || haskey(kwargs,:x0)
        return MichelsenTPFlash(;kwargs...)
    elseif haskey(kwargs,:numphases)
        return DETPFlash(;kwargs...)
    else
        return DETPFlash(;kwargs...)
    end
end

export tp_flash
