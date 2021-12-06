"""
    TPFlashMethod

Abstract type for `tp_flash` routines. it requires defining `numphases(method)` and `tp_flash_impl(model,p,T,n,method)`.

"""
abstract type TPFlashMethod end

"""
tp_flash(model, p, T, n, method::TPFlashMethod =DETPFlash())

Routine to solve non-reactive multicomponent flash problem.
The default method uses Global Optimization. see [DETPFlash](@ref)

Inputs: 
 - T, Temperature
 - p, Pressure
 - n, vector of number of moles of each species

Outputs - Tuple containing: 
 - xᵢⱼ, Array of mole fractions of species j in phase i
 - nᵢⱼ, Array of mole numbers of species j in phase i, [mol]
 - G, Gibbs Free Energy of Equilibrium Mixture [J]
"""
function tp_flash end

"""
    numphases(method::TPFlashMethod)

return the number of phases supported by the TP flash method. by default its set to 2.
it the method allows it, you can set the number of phases by doing `method(;numphases = n)`.
"""
numphases(method::TPFlashMethod) = 2

#default
include("tp_flash/DifferentialEvolutiontp_flash.jl")
include("tp_flash/RachfordRicetp_flash.jl")

function tp_flash(model::EoSModel, p, T, n,method::TPFlashMethod = DETPFlash())
    numspecies = length(model)
    if numspecies != length(n)
        error("There are ", numspecies,
            " species in the model, but the number of mole numbers specified is ", 
            length(n))
    end
    if numphases(method) == 1
        return (n, n / sum(n), gibbs_free_energy(model, p, T, n))
    end
    return tp_flash_impl(model::EoSModel,p,T,n,method)
end

export tp_flash
