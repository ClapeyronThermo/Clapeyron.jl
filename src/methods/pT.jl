const VT_STRING = 
"""

For Helmholtz-based models, it calls [`Clapeyron.volume`](@ref) to obtain `V` and evaluate the property in a volume-temperature (VT) basis.
The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
If you already have a volume, the [`VT0`](@ref) module is available to evaluate this property directly in the VT basis, bypassing the volume iterative calculation.

Gibbs-based models are instead evaluated directly in the pressure-temperature basis.

"""

const IDEALMODEL_REQUIRED = 
"""
!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel = ReidIdeal)`).
"""

const SINGLE_PHASE_PROP = 
"""
!!! note "single phase property"
    This property is not defined for more than one phase. Calling this property with two-phase states will result in an error.
"""

"""
    volume(model::EoSModel, p, T, z=SA[1.0]; phase=:unknown, threaded=true, vol0=nothing)
    volume(result::FlashResult)
    volume(model, result::FlashResult)
    volume(result::FlashResult, phase_index::Int)
    volume(model, result::FlashResult, phase_index::Int)


Calculates the volume `(m³)` of the compound modelled by `model` at a certain pressure `p`, temperature `T` and moles `z`.
`phase` is a Symbol that determines the initial volume root to look for:
- If `phase =:unknown` (Default), it will return the physically correct volume root with the least Gibbs energy.
- If `phase =:liquid`, it will return the volume of the phase using a liquid initial point.
- If `phase =:vapor`, it will return the volume of the phase using a gas initial point.
- If `phase =:solid`, it will return the volume of the phase using a solid initial point (only supported for EoS that support a solid phase).
- If `phase =:stable`, it will return the physically correct volume root with the least Gibbs energy, and perform a stability test on the result.

All volume calculations are checked for mechanical stability, that is: `dP/dV <= 0`.

The calculation of both volume roots can be calculated in serial (`threaded=false`) or in parallel (`threaded=true`).

An initial estimate of the volume `vol0` can be optionally be provided.

`volume(result::FlashResult)` will return the volume of the aggregate of phases stored in the `FlashResult` whereas `volume(result::FlashResult,phase_index)` will return the volume of the ith phase. 
Because molar volumes are directly stored in the `FlashResult` struct, `volume(model,result)` will just call `volume(result)` instead.
Similarly, `volume(model,result::FlashResult,i)` will just call `volume(result,i)`.

!!! tip
    The volume computation may fail and return `NaN` because the default initial point is too far from the actual volume.
    Providing a value for `vol0` may help in these situations.
    Such a starting point can be found from physical knowledge, or by computing the volume using a different model for example.

!!! warning "Stability checks"
    The stability check is disabled by default. That means that the volume obtained just follows the relation `p = pressure(model,V,T,z)`.
    For single component models, this is alright, but phase splits (with different compositions that the input) can and will occur, meaning that
    the volume solution does not correspond to an existing phase.
    For unknown multicomponent mixtures, it is recommended to use a phase equilibrium procedure (like `tp_flash`) to obtain a list of valid compositions, and then perform a volume calculation over those compositions.
    You can also pass `phase=:stable` to perform the stability test inside the volume solver. Finally, you can perform the stability test after the volume solver:
    ```julia
    v = volume(model,p,T,z)
    isstable(model,v,T,z)
    ```
"""
function volume end

"""
    molar_density(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    molar_density(result::FlashResult)
    molar_density(model, result::FlashResult)
    molar_density(result::FlashResult, phase_index::Int)
    molar_density(model, result::FlashResult, phase_index::Int)

Default units: `[mol·m⁻³]`

Calculates the molar density, defined as:

```julia
ρₙ = ∑nᵢ/V
```

`molar_density(model,result::FlashResult)` will return the molar density of the aggregate of phases stored in the `FlashResult` whereas `molar_density(model,result::FlashResult,i::Int)` will return the molar density of the ith phase. 
Because molar volumes are directly stored in the `FlashResult` struct, `molar_density(model,result)` will just call `molar_density(result)` instead.
Similarly, `molar_density(model,result::FlashResult,i)` will just call `molar_density(result,i)`.

$VT_STRING
"""
function molar_density(model::EoSModel,p,T,z=SA[1.0];phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_molar_density(model,V,T,z)
end

"""
    mass_density(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true)
    mass_density(model, result::FlashResult)
    mass_density(model, result::FlashResult, phase_index::Int)

Default units: `[kg·m⁻³]`

Calculates the mass density, defined as:

```julia
ρₙ = Mr/V
```
Where `Mr` is the molecular weight of the model at the input composition.

`mass_density(model,result::FlashResult)` will return the mass density of the aggregate of phases stored in the `FlashResult` whereas `mass_density(model,result::FlashResult,i::Int)` will return the mass density of the ith phase. 


$VT_STRING
"""
function mass_density(model::EoSModel, p, T, z=SA[1.0]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_mass_density(model,V,T,z)
end

"""
    compressibility_factor(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    compressibility_factor(result::FlashResult)
    compressibility_factor(model, result::FlashResult)
    compressibility_factor(result::FlashResult, phase_index::Int)
    compressibility_factor(model, result::FlashResult, phase_index::Int)

Calculates the compressibility factor `Z`, defined as:

```julia
Z = p*V(p)/R*T
```

`compressibility_factor(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`compressibility_factor(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.
Because molar volumes and pressures are directly stored in the `FlashResult` struct, `compressibility_factor(model,result)` will just call `compressibility_factor(result)` instead.
Similarly, `compressibility_factor(model,result::FlashResult,i)` will just call `compressibility_factor(result,i)`.

$VT_STRING
$SINGLE_PHASE_PROP
"""
function compressibility_factor(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    #this property only depends on the implementation of volume_impl.
    V = volume(model,p,T,z;phase,threaded,vol0)
    return p*V/(sum(z)*Rgas(model)*T)
end

function PT_property(model,p,T,z,phase,threaded,vol0,f::F,vol::VV) where {F,VV}
    
    if f == pressure
        return p
    elseif f == temperature
        return T
    end

    if z isa Number
        return PT_property(model,p,T,SA[z],phase,threaded,vol0,f)
    end
    if isnothing(vol)
        V = volume(model, p, T, z; phase, threaded, vol0)
    else
        V = one(Base.promote_eltype(model,p,T,z))*vol
    end
    if VT_use_p(f)
        return f(model,V,T,z,p)
    else
        return f(model,V,T,z)
    end
end

PT_property(model,p,T,z,phase,threaded,vol0,f::F) where {F} = PT_property(model,p,T,z,phase,threaded,vol0,f,nothing)
PT_property(model,p,T,z,phase,vol,f::F) where {F} = PT_property(model,p,T,z,phase,false,nothing,f,vol)
"""
    entropy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    entropy(model, result::FlashResult)
    entropy(model, result::FlashResult, phase_index::Int)

Default units: `[J·K⁻¹]`

Calculates entropy, defined as:

```julia
S = -∂A/∂T
```

`entropy(model,result::FlashResult)` will return the entropy of the aggregate of phases stored in the `FlashResult` whereas `entropy(model,result::FlashResult,i::Int)` will return the entropy of the ith phase. 

$VT_STRING
"""
function entropy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    return PT_property(model,p,T,z,phase,threaded,vol0,VT_entropy)
end

"""
    mass_entropy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    mass_entropy(model, result::FlashResult)
    mass_entropy(model, result::FlashResult, phase_index::Int)

Default units: `[J·kg⁻¹·K⁻¹]`

Calculates entropy per unit of mass, defined as:

```julia
S = -∂A/∂T/Mr
```
Where `Mr` is the molecular weight of the model at the input composition.

`mass_entropy(model,result::FlashResult)` will return the mass entropy of the aggregate of phases stored in the `FlashResult` whereas `mass_entropy(model,result::FlashResult,i::Int)` will return the mass entropy of the ith phase. 

$VT_STRING
"""
function mass_entropy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    return PT_property(model,p,T,z,phase,threaded,vol0,VT_mass_entropy)
end

"""
    entropy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    entropy_res(model, result::FlashResult)
    entropy_res(model, result::FlashResult, phase_index::Int)

Default units: `[J·K⁻¹]`

Calculates residual entropy, defined as:

```julia
S = -∂Ares/∂T
```

`entropy_res(model,result::FlashResult)` will return the residual entropy of the aggregate of phases stored in the `FlashResult` whereas `entropy_res(model,result::FlashResult,i::Int)` will return the residual entropy of the ith phase. 

$VT_STRING
"""
function entropy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    return PT_property(model,p,T,z,phase,threaded,vol0,VT_entropy_res)
end

"""
    chemical_potential(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    chemical_potential(model, result::FlashResult)
    chemical_potential(model, result::FlashResult, phase_index::Int)

Default units: `[J·mol⁻¹]`

Calculates the chemical potential, defined as:

```julia
μᵢ = ∂A/∂nᵢ
```

`chemical_potential(model,result::FlashResult)` will return the chemical potential of the aggregate of phases stored in the `FlashResult` whereas `chemical_potential(model,result::FlashResult,i::Int)` will return the chemical potential of the ith phase.

$VT_STRING
"""
function chemical_potential(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    μ = chemical_potential_impl(model,p,T,z,phase,threaded,vol0)
end

function chemical_potential_impl(model,p,T,z,phase,threaded,vol0)
    return PT_property(model,p,T,z,phase,threaded,vol0,VT_chemical_potential)
end

"""
    chemical_potential_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    chemical_potential_res(model, result::FlashResult)
    chemical_potential_res(model, result::FlashResult, phase_index::Int)

Default units: `[J·mol⁻¹]`

Calculates the residual chemical potential, defined as:

```julia
μresᵢ = ∂Ares/∂nᵢ
```

`chemical_potential_res(model,result::FlashResult)` will return the residual chemical potential of the aggregate of phases stored in the `FlashResult` whereas `chemical_potential_res(model,result::FlashResult,i::Int)` will return the residual chemical potential of the ith phase.

$VT_STRING
"""
function chemical_potential_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_chemical_potential_res)
end

"""
    internal_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    internal_energy(model, result::FlashResult)
    internal_energy(model, result::FlashResult, phase_index::Int)

Default units: `[J]`

Calculates the internal energy, defined as:

```julia
U = A - T * ∂A/∂T
```

`internal_energy(model,result::FlashResult)` will return the internal energy of the aggregate of phases stored in the `FlashResult` whereas `internal_energy(model,result::FlashResult,i::Int)` will return the internal energy of the ith phase.

$VT_STRING
"""
function internal_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_internal_energy)
end

"""
    mass_internal_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    mass_internal_energy(model, result::FlashResult)
    mass_internal_energy(model, result::FlashResult, phase_index::Int)

Default units: `[J·kg⁻¹]`

Calculates the internal energy per unit of mass, defined as:

```julia
U = (A - T * ∂A/∂T)/Mr
```
Where `Mr` is the molecular weight of the model at the input composition.

`mass_internal_energy(model,result::FlashResult)` will return the mass internal energy of the aggregate of phases stored in the `FlashResult` whereas `mass_internal_energy(model,result::FlashResult,i::Int)` will return the mass internal energy of the ith phase.

$VT_STRING
"""
function mass_internal_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_mass_internal_energy)
end

"""
    internal_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    internal_energy_res(model, result::FlashResult)
    internal_energy_res(model, result::FlashResult, phase_index::Int)

Default units: `[J]`

Calculates the residual internal energy, defined as:

```julia
U = Ar - T * ∂Ar/∂T
```

`internal_energy_res(model,result::FlashResult)` will return the residual internal energy of the aggregate of phases stored in the `FlashResult` whereas `internal_energy_res(model,result::FlashResult,i::Int)` will return the residual internal energy of the ith phase.

$VT_STRING
"""
function internal_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_internal_energy_res)
end

"""
    enthalpy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    enthalpy(model, result::FlashResult)
    enthalpy(model, result::FlashResult, phase_index::Int)

Default units: `[J]`

Calculates the enthalpy, defined as:

```julia
H = A - T * ∂A/∂T - V * ∂A/∂V
```

`enthalpy(model,result::FlashResult)` will return the enthalpy of the aggregate of phases stored in the `FlashResult` whereas `enthalpy(model,result::FlashResult,i::Int)` will return the enthalpy of the ith phase.

$VT_STRING
"""
function enthalpy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_enthalpy)
end

"""
    mass_enthalpy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    mass_enthalpy(model, result::FlashResult)
    mass_enthalpy(model, result::FlashResult, phase_index::Int)

Default units: `[J·kg⁻¹]`

Calculates the enthalpy per unit of mass, defined as:

```julia
H = (A - T * ∂A/∂T - V * ∂A/∂V)/Mr
```
Where `Mr` is the molecular weight of the model at the input composition.

`mass_enthalpy(model,result::FlashResult)` will return the mass enthalpy of the aggregate of phases stored in the `FlashResult` whereas `mass_enthalpy(model,result::FlashResult,i::Int)` will return the mass enthalpy of the ith phase.

$VT_STRING
"""
function mass_enthalpy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_mass_enthalpy)
end

"""
    enthalpy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    enthalpy_res(model, result::FlashResult)
    enthalpy_res(model, result::FlashResult, phase_index::Int)

Default units: `[J]`

Calculates the residual enthalpy, defined as:

```julia
H = Ar - T * ∂Ar/∂T - V * ∂Ar/∂V
```

`enthalpy_res(model,result::FlashResult)` will return the residual enthalpy of the aggregate of phases stored in the `FlashResult` whereas `enthalpy_res(model,result::FlashResult,i::Int)` will return the residual enthalpy of the ith phase.

$VT_STRING
"""
function enthalpy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_enthalpy_res)
end

"""
    gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    gibbs_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    gibbs_free_energy(model, result::FlashResult)
    gibbs_free_energy(model, result::FlashResult, phase_index::Int)

Default units: `[J]`

Calculates the Gibbs energy, defined as:

```julia
G = A + p*V
```

`gibbs_free_energy(model,result::FlashResult)` will return the Gibbs energy of the aggregate of phases stored in the `FlashResult` whereas `gibbs_free_energy(model,result::FlashResult,i::Int)` will return the Gibbs energy of the ith phase.

$VT_STRING
"""
function gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_gibbs_free_energy)
end

"""
    mass_gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    mass_gibbs_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    mass_gibbs_free_energy(model, result::FlashResult)
    mass_gibbs_free_energy(model, result::FlashResult, phase_index::Int)

Default units: `[J·kg⁻¹]`

Calculates the Gibbs energy per unit of mass, defined as:

```julia
G = (A + p*V)/Mr
```
Where `Mr` is the molecular weight of the model at the input composition.

`mass_gibbs_free_energy(model,result::FlashResult)` will return the mass Gibbs energy of the aggregate of phases stored in the `FlashResult` whereas `mass_gibbs_free_energy(model,result::FlashResult,i::Int)` will return the mass Gibbs energy of the ith phase.

$VT_STRING
"""
function mass_gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_mass_gibbs_free_energy)
end

"""
    gibbs_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    gibbs_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    gibbs_free_energy_res(model, result::FlashResult)
    gibbs_free_energy_res(model, result::FlashResult, phase_index::Int)

Default units: `[J]`

Calculates the residual Gibbs energy, defined as:

```julia
G = Ar - V*∂Ar/∂V
```

`gibbs_free_energy_res(model,result::FlashResult)` will return the residual Gibbs energy of the aggregate of phases stored in the `FlashResult` whereas `gibbs_free_energy_res(model,result::FlashResult,i::Int)` will return the residual Gibbs energy of the ith phase.

$VT_STRING
"""
function gibbs_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_gibbs_free_energy_res)
end

"""
    helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    helmholtz_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    helmholtz_free_energy(model, result::FlashResult)
    helmholtz_free_energy(model, result::FlashResult, phase_index::Int)

Default units: `[J]`

Calculates the Helmholtz energy, defined as:

```julia
A = eos(model,V(p),T,z)
```

`helmholtz_free_energy(model,result::FlashResult)` will return the Helmholtz energy of the aggregate of phases stored in the `FlashResult` whereas `helmholtz_free_energy(model,result::FlashResult,i::Int)` will return the Helmholtz energy of the ith phase.

$VT_STRING
"""
function helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_helmholtz_free_energy)
end

"""
    mass_helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    mass_helmholtz_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    mass_helmholtz_free_energy(model, result::FlashResult)
    mass_helmholtz_free_energy(model, result::FlashResult, phase_index::Int)

Default units: `[J·kg⁻¹]`

Calculates the Helmholtz energy per unit of mass, defined as:

```julia
A = eos(model,V(p),T,z)/Mr
```
Where `Mr` is the molecular weight of the model at the input composition.

`mass_helmholtz_free_energy(model,result::FlashResult)` will return the mass Helmholtz energy of the aggregate of phases stored in the `FlashResult` whereas `mass_helmholtz_free_energy(model,result::FlashResult,i::Int)` will return the mass Helmholtz energy of the ith phase.

$VT_STRING
"""
function mass_helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_mass_helmholtz_free_energy)
end

"""
    helmholtz_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    helmholtz_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    helmholtz_free_energy_res(model, result::FlashResult)
    helmholtz_free_energy_res(model, result::FlashResult, phase_index::Int)

Default units: `[J]`

Calculates the residual Helmholtz energy, defined as:

```julia
A = eos_res(model,V(p),T,z)
```

`helmholtz_free_energy_res(model,result::FlashResult)` will return the residual Helmholtz energy of the aggregate of phases stored in the `FlashResult` whereas `helmholtz_free_energy_res(model,result::FlashResult,i::Int)` will return the residual Helmholtz energy of the ith phase.

$VT_STRING
"""
function helmholtz_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_helmholtz_free_energy_res)
end

const helmholtz_energy = helmholtz_free_energy
const helmholtz_energy_res = helmholtz_free_energy_res 
const gibbs_energy = gibbs_free_energy
const gibbs_energy_res = gibbs_free_energy_res
const mass_helmholtz_energy = mass_helmholtz_free_energy
const mass_gibbs_energy = mass_gibbs_free_energy
"""
    isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    isochoric_heat_capacity(model, result::FlashResult)
    isochoric_heat_capacity(model, result::FlashResult, phase_index::Int)

Default units: `[J·K⁻¹]`

Calculates the isochoric heat capacity, defined as:

```julia
Cv = -T * ∂²A/∂T²
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isochoric_heat_capacity(model,V,T,z)`.

`isochoric_heat_capacity(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`isochoric_heat_capacity(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING

"""
function isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isochoric_heat_capacity)
end

"""
    mass_isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    mass_isochoric_heat_capacity(model, result::FlashResult)
    mass_isochoric_heat_capacity(model, result::FlashResult, phase_index::Int)

Default units: `[J·kg⁻¹·K⁻¹]`

Calculates the isochoric heat capacity per unit of mass, defined as:

```julia
Cv = -T * ∂²A/∂T² / Mr
```
Where `Mr` is the molecular weight of the model at the input composition.

`mass_isochoric_heat_capacity(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`mass_isochoric_heat_capacity(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
$IDEALMODEL_REQUIRED
$SINGLE_PHASE_PROP
"""
function mass_isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_mass_isochoric_heat_capacity)
end

"""
    isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    isobaric_heat_capacity(model, result::FlashResult)
    isobaric_heat_capacity(model, result::FlashResult, phase_index::Int)

Default units: `[J·K⁻¹]`

Calculates the isobaric heat capacity, defined as:

```julia
Cp = -T*(∂²A/∂T² - (∂²A/∂V∂T)^2 / ∂²A/∂V²)
```

`isobaric_heat_capacity(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`isobaric_heat_capacity(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
$IDEALMODEL_REQUIRED
$SINGLE_PHASE_PROP
"""
function isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isobaric_heat_capacity)
end

"""
    mass_isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    mass_isobaric_heat_capacity(model, result::FlashResult)
    mass_isobaric_heat_capacity(model, result::FlashResult, phase_index::Int)

Default units: `[J·kg⁻¹·K⁻¹]`

Calculates the isobaric heat capacity per unit of mass, defined as:

```julia
Cp = (-T*(∂²A/∂T² - (∂²A/∂V∂T)^2 / ∂²A/∂V²))/Mr
```
Where `Mr` is the molecular weight of the model at the input composition.

`mass_isobaric_heat_capacity(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`mass_isobaric_heat_capacity(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
$IDEALMODEL_REQUIRED
$SINGLE_PHASE_PROP
"""
function mass_isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_mass_isobaric_heat_capacity)
end

"""
    adiabatic_index(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    adiabatic_index(model, result::FlashResult)
    adiabatic_index(model, result::FlashResult, phase_index::Int)

Calculates the adiabatic index, defined as:

```julia
γ = Cp/Cv
```

`adiabatic_index(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`adiabatic_index(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
$IDEALMODEL_REQUIRED
$SINGLE_PHASE_PROP
"""
function adiabatic_index(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_adiabatic_index)
end

"""
    isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    isothermal_compressibility(model, result::FlashResult)
    isothermal_compressibility(model, result::FlashResult, phase_index::Int)

Default units: `[Pa⁻¹]`

Calculates the isothermal compressibility, defined as:

```julia
κₜ = -(V*∂p/∂V)⁻¹
```

`isothermal_compressibility(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`isothermal_compressibility(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
$SINGLE_PHASE_PROP
"""
function isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isothermal_compressibility)
end

"""
    isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    isentropic_compressibility(model, result::FlashResult)
    isentropic_compressibility(model, result::FlashResult, phase_index::Int)

Default units: `[Pa⁻¹]`

Calculates the isentropic compressibility, defined as:

```julia
κₛ = (V*( ∂²A/∂V² - ∂²A/∂V∂T^2 / ∂²A/∂T² ))⁻¹
```

`isentropic_compressibility(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`isentropic_compressibility(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
$IDEALMODEL_REQUIRED
$SINGLE_PHASE_PROP
"""
function isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isentropic_compressibility)
end

"""
    speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    speed_of_sound(model, result::FlashResult)
    speed_of_sound(model, result::FlashResult, phase_index::Int)

Default units: `[m·s⁻¹]`

Calculates the speed of sound, defined as:

```julia
c = V * √(∂²A/∂V² - ∂²A/∂V∂T^2 / ∂²A/∂T²)/Mr)
```
Where `Mr` is the molecular weight of the model at the input composition.

`speed_of_sound(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`speed_of_sound(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
$IDEALMODEL_REQUIRED
$SINGLE_PHASE_PROP
"""
function speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_speed_of_sound)
end

"""
    isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    isobaric_expansivity(model, result::FlashResult)
    isobaric_expansivity(model, result::FlashResult, phase_index::Int)

Default units: `[K⁻¹]`

Calculates the isobaric expansivity, defined as:

```julia
α = -∂²A/∂V∂T / (V*∂²A/∂V²)
```

`isobaric_expansivity(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`isobaric_expansivity(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
$SINGLE_PHASE_PROP
"""
function isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isobaric_expansivity)
end

"""
    joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    joule_thomson_coefficient(model, result::FlashResult)
    joule_thomson_coefficient(model, result::FlashResult, phase_index::Int)

Default units: `[K·Pa⁻¹]`

Calculates the Joule–Thomson coefficient, defined as:

```julia
μⱼₜ = -(∂²A/∂V∂T - ∂²A/∂V² * ((T*∂²A/∂T² + V*∂²A/∂V∂T) / (T*∂²A/∂V∂T + V*∂²A/∂V²)))⁻¹
```

`joule_thomson_coefficient(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`joule_thomson_coefficient(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
$IDEALMODEL_REQUIRED
$SINGLE_PHASE_PROP
"""
function joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_joule_thomson_coefficient)
end

"""
    identify_phase(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)::Symbol
    identify_phase(result::FlashResult, i::Int)
    identify_phase(model::EoSModel, result::FlashResult, i::Int)

Returns the phase of a fluid at the conditions specified by `V`, `T` and `z`.
Uses the phase identification parameter criteria from `Clapeyron.pip`.

Returns `:liquid` if the phase is liquid (or liquid-like), `:vapour` if the phase is vapour (or vapour-like), and `:unknown` if the calculation of the phase identification parameter failed.

`identify_phase(model,result::FlashResult,i::Int)` will return the phase type of the ith phase stored in the result, if available. 
`identify_phase(model,result,i)` will try to get the stored result first, and fallback to `identify_phase(model,p,T,xi)` if there is no phase stored.

$VT_STRING
"""
function identify_phase(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, vol = NaN)
    if isnan(vol) || isnothing(vol)
        V = volume(model, p, T, z; phase, threaded, vol0)
    else
        V = vol*oneunit(Base.promote_eltype(model,p,T,z))
    end
    return VT_identify_phase(model,V,T,z)
end


"""
    fundamental_derivative_of_gas_dynamics(model::EoSModel, p, T, z=SA[1.]; phase=:gas, threaded=true, vol0=nothing)::Symbol
    fundamental_derivative_of_gas_dynamics(model, result::FlashResult)
    fundamental_derivative_of_gas_dynamics(model, result::FlashResult, phase_index::Int)

Calculates the fundamental derivative of gas dynamics.

`fundamental_derivative_of_gas_dynamics(model,result::FlashResult)` will calculate the property only if there is a single phase in the result, and error otherwise.
`fundamental_derivative_of_gas_dynamics(model,result::FlashResult, i::Int)` will calculate the property at the ith phase of a `FlashResult`.

$VT_STRING
"""
function fundamental_derivative_of_gas_dynamics(model::EoSModel, p, T, z=SA[1.]; phase=:gas, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_fundamental_derivative_of_gas_dynamics)
end

"""
    fugacity_coefficient(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Calculates the fugacity coefficient φᵢ, defined as:

```julia
log(φᵢ) = μresᵢ/RT - log(Z)
```
Where `μresᵢ` is the vector of residual chemical potentials and `Z` is the compressibility factor.

$VT_STRING
"""
function fugacity_coefficient(model::EoSModel,p,T,z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_fugacity_coefficient)
end

function fugacity_coefficient!(φ,model::EoSModel,p,T,z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    VT_fugacity_coefficient!(φ,model,V,T,z,p)
end

"""
    activity_coefficient(model::EoSModel,p,T,z=SA[1.0];reference = :pure, phase=:unknown, threaded=true, vol0=nothing)

Calculates the activity, defined as:

```julia
log(γ*z) = (μ_mixt - μ_ref) / R̄ / T
```
where `μ_mixt` is the chemical potential of the mixture and `μ_ref` is the reference chemical potential for the model at `p`,`T` conditions, calculated via [`Clapeyron.reference_chemical_potential`](@ref).
If the `μ_ref` keyword argument is not provided, the `reference` keyword is used to specify the reference chemical potential.

$VT_STRING

"""
function activity_coefficient(model::EoSModel,p,T,z = SA[1.0];
                            μ_ref = nothing,
                            reference = :pure,
                            phase=:unknown,
                            threaded=true,
                            vol0=nothing)
    reference = Symbol(reference)
    γmodel = __γ_unwrap(model)
    if γmodel isa ActivityModel
        return activity_coefficient(γmodel,p,T,z)
    end
    if μ_ref == nothing
        return activity_coefficient_impl(model,p,T,z,reference_chemical_potential(model,p,T,reference;phase,threaded),reference,phase,threaded,vol0)
    else
        return activity_coefficient_impl(model,p,T,z,μ_ref,reference,phase,threaded,vol0)
    end
end

function activity_coefficient_impl(model,p,T,z,μ_ref,reference,phase,threaded,vol0)
    RT = Rgas(model)*T
    μ_mixt = chemical_potential(model, p, T, z; phase, threaded, vol0)
    return sum(z) .* exp.((μ_mixt .- μ_ref) ./ RT) ./z
end

"""
    activity(model::EoSModel,p,T,z=SA[1.0];reference = :pure, phase=:unknown, threaded=true, vol0=nothing)

Calculates the activity, defined as:

```julia
log(a) = (μ_mixt - μ_ref) / R̄ / T
```
where `μ_mixt` is the chemical potential of the mixture and `μ_ref` is the reference chemical potential for the model at `p`,`T` conditions, calculated via [`Clapeyron.reference_chemical_potential`](@ref).
If the `μ_ref` keyword argument is not provided, the `reference` keyword is used to specify the reference chemical potential.

$VT_STRING
"""
function activity(model::EoSModel,p,T,z;
                μ_ref = nothing,
                reference = :pure,
                phase=:unknown,
                threaded=true,
                vol0=nothing)
    reference = Symbol(reference)
    if model isa ActivityModel
        return activity(model,p,T,z)
    end
    if μ_ref == nothing
        return activity_impl(__γ_unwrap(model),p,T,z,reference_chemical_potential(model,p,T,reference;phase,threaded),reference,phase,threaded,vol0)
    else
        return activity_impl(__γ_unwrap(model),p,T,z,μ_ref,reference,phase,threaded,vol0)
    end
end

function activity_impl(model,p,T,z,μ_ref,reference,phase,threaded,vol0)
    R̄ = Rgas(model)
    μ_mixt = chemical_potential(model, p, T, z; phase, threaded, vol0)
    return exp.((μ_mixt .- μ_ref) ./ R̄ ./ T)
end

function find_hydronium_index(model)
    idx = findfirst(isequal("hydronium"),component_list(model))
    idx == nothing && return 0
    return idx
end

function find_hydroxide_index(model)
    idx = findfirst(isequal("hydroxide"),component_list(model))
    idx == nothing && return 0
    return idx
end

function find_water_indx(model)
    idx = findfirst(isequal("water"),component_list(model))
    idx == nothing && return 0
    return idx
end

"""
    aqueous_activity(model::EoSModel,p,T,z=SA[1.0]; phase=:unknown, threaded=true, vol0=nothing)

Calculates the activity with the reference being infinite dilution in water, defined as:

```julia
log(a) = (μ_mixt - μ_inf) / R̄ / T
```
where `μ_mixt` is the chemical potential of the mixture and `μ_inf` is the chemical potential of the components at infinite dilution in water.

$VT_STRING
"""
function aqueous_activity(model::EoSModel,p,T,z=SA[1.];                            
                        μ_ref = nothing,
                        reference = :aqueous,
                        phase=:unknown,
                        threaded=true,
                        vol0=nothing)
    return activity(model,p,T,z;reference=reference,μ_ref = μ_ref,phase=phase,threaded=threaded,vol0=vol0)
end


"""
    reference_chemical_potential_type(model)::Symbol

Returns a symbol with the type of reference chemical potential used by the input `model`.

"""
reference_chemical_potential_type(model) = :pure

"""
    reference_chemical_potential(model::EoSModel,p,T,reference; phase=:unknown, threaded=true, vol0=nothing)

Returns a reference chemical potential. Used in calculation of `activity` and activity_coefficient. There are two available references:
- `:pure`: the reference potential is a pure component at specified `T`, `p` and `phase`
- `:aqueous`: the chemical potential of the pure components at specified `T`, `p` and `phase`
- `:sat_pure_T`:  the reference potential is the pure saturated liquid phase at specified `T`.
- `:zero`: the reference potential is equal to zero for all components (used for `ActivityModel`)

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function reference_chemical_potential(model::EoSModel,p,T,reference = reference_chemical_potential_type(model); phase=:unknown, threaded=true, vol0=nothing)
    reference = Symbol(reference)
    if reference == :pure
        pure = split_pure_model(model)
        return gibbs_free_energy.(pure, p, T; phase, threaded)
    elseif reference == :aqueous
        idx_w = find_water_indx(model)
        if idx_w == 0
            throw(ArgumentError("There is no water in $(model)."))
        end
        zref = ones(length(model))
        zref[1:length(model) .!= idx_w] .*= 0.01801528
        zref ./= sum(zref)
        return chemical_potential(model, p, T, zref; phase, threaded, vol0)
    elseif reference == :sat_pure_T
        pure = split_pure_model(model)
        sat = saturation_pressure.(pure,T)
        vl_pure = getindex.(sat,2)
        return VT_gibbs_free_energy.(pure, vl_pure, T)
    elseif reference == :zero
        _0 = Base.promote_eltype(model,p,T)
        return fill(_0,length(model))
    else
        throw(ArgumentError("reference must be one of :pure, :aqueous, :sat_pure_T, :zero"))
    end
end

"""
    inversion_temperature(model::EoSModel, p, z=SA[1.0]; phase=:unknown, threaded=true, vol0=nothing)

Calculates the inversion temperature `T_inv`, defined as the temperature where the Joule-Thomson coefficient becomes zero, i.e.

```julia
μⱼₜ(T) = 0
```
The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

See also [`joule_thomson_coefficient`](@ref).
"""
function inversion_temperature(model::EoSModel, p, z=SA[1.0]; phase=:unknown, threaded=true, vol0=nothing)
    T0 = 6.75*T_scale(model,z)
    μⱼₜ(T) = joule_thomson_coefficient(model, p, T, z; phase, threaded, vol0)
    return Roots.find_zero(μⱼₜ,T0)
end


"""
    mixing(model::EoSModel, p, T, z=SA[1.], property; phase=:unknown, threaded=true, vol0=nothing)

Calculates the mixing function for a specified property as:

```julia
f_mix = f(p,T,z) - ∑zᵢ*f_pureᵢ(p,T)
```
The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function mixing(model::EoSModel, p, T, z, property::ℜ; phase=:unknown, threaded=true, vol0=nothing, output=nothing) where {ℜ}
    pure = split_pure_model(model)
    TT = typeof(p+T+first(z))
    mix_prop  = property(model, p, T, z; phase, threaded, vol0)
    for i in 1:length(z)
        mix_prop -= z[i]*property(pure[i], p, T; phase, threaded)
    end
    return mix_prop::TT
end

"""
    excess(model::EoSModel, p, T, z, property; phase=:unknown, threaded=true, vol0=nothing)

Returns the excess value of a bulk property relative to its ideal mixing value.

By default this delegates to [`mixing`](@ref). For some properties (e.g.
`entropy` and `gibbs_free_energy`) specialized implementations are provided to
use residual contributions.
"""
function excess(model::EoSModel, p, T, z, property; phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    mixing(model, p, T, z, property; phase, threaded, vol0, output)
end

function excess(model::EoSModel, p, T, z, ::typeof(entropy); phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    TT = typeof(p+T+first(z))
    pure = split_pure_model(model)
    s_mix = entropy_res(model, p, T, z; phase, threaded, vol0)
    for i in 1:length(z)
        s_mix -= z[i]*entropy_res(pure[i], p, T; phase, threaded)
    end
    #s_pure = entropy_res.(pure,p,T)
    return s_mix::TT
end

function excess(model::EoSModel, p, T, z, ::typeof(gibbs_free_energy); phase=:unknown, threaded=true, vol0=nothing, output=nothing)
    TT = typeof(p+T+first(z))
    pure = split_pure_model(model)
    g_mix = gibbs_free_energy(model, p, T, z; phase, threaded, vol0)
    log∑z = log(sum(z))
    R̄ = Rgas(model)
    for i in 1:length(z)
        lnxi = R̄*T*(log(z[i]) - log∑z)
        g_mix -= z[i]*(gibbs_free_energy(pure[i], p, T; phase, threaded) + lnxi)
    end

    return g_mix::TT
end


"""
    gibbs_solvation(model::EoSModel, T; threaded=true, vol0=(nothing,nothing))

Calculates the solvation Gibbs energy as:

```julia
g_solv = -R̄*T*log(K)
```
Where the first component is the solvent and second is the solute.
"""
function gibbs_solvation(model::EoSModel, T; threaded=true, vol0=(nothing,nothing))
    binary_component_check(gibbs_solvation, model)
    pure = split_pure_model(model)
    z = [1.0,1e-30]

    p,v_l,v_v = saturation_pressure(pure[1],T)

    φ_l = fugacity_coefficient(model, p, T, z; phase=:l, threaded, vol0=vol0[1])
    φ_v = fugacity_coefficient(model, p, T, z; phase=:v, threaded, vol0=vol0[2])

    K = φ_v[2]*v_v/φ_l[2]/v_l
    R̄ = Rgas(model)
    return -R̄*T*log(K)
end

"""
    partial_property(model::EoSModel, p, T, z, property::X; phase=:unknown, threaded=true, vol0=nothing) where {X} is any extensive property.

Calculates the partial molar property of a mixture at specified temperature, pressure, mol amounts, and extensive property of interest.
The equality `sum(z .* partial_property(model,p,T,z,property) - property(model,p,T,z))` should hold.
    
The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function partial_property(model::EoSModel, p, T, z, property::ℜ; phase=:unknown, threaded=true, vol0=nothing, output=nothing) where {ℜ}
    V = volume(model, p, T, z; phase, threaded, vol0)
    return _partial_property(model,V,T,z,PT_to_VT(property))
end

function _partial_property(model::EoSModel, V, T, z::AbstractVector, ::typeof(volume))
    _,∂p∂V = p∂p∂V(model,V,T,z)
    ∂p∂nᵢ = VT_molar_gradient(model, V, T, z, pressure)
    return -∂p∂nᵢ ./ ∂p∂V
end

function _partial_property(model::EoSModel, V, T, z::AbstractVector, ::typeof(VT_gibbs_free_energy))
    return VT_molar_gradient(model,V,T,z,eos)
end

function _partial_property(model::EoSModel, V, T, z::AbstractVector, VT_prop::F) where F
    ∂x∂nᵢ = VT_molar_gradient(model,V,T,z,VT_prop)
    #triple product rule:
    #∂x∂nᵢ|p = ∂x∂nᵢ|V - ∂x∂V * ∂p∂nᵢ|V * ∂p∂V^-1
    ∂p∂nᵢ = VT_molar_gradient(model,V,T,z,pressure)
    xv = @deferred_V(VT_prop,partial_property)
    ∂x∂V = Solvers.derivative(xv,V)
    _,∂p∂V = p∂p∂V(model,V,T,z)
    return ∂x∂nᵢ .- ∂x∂V .* ∂p∂nᵢ ./ ∂p∂V
end

"""
    thermodynamic_factor(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Calculates the thermodynamic factor matrix Γᵢⱼ (size: N-1 × N-1) defined as:

```julia
Γᵢⱼ = δᵢⱼ + xᵢ ∂lnγᵢ/∂xⱼ
```
"""
function thermodynamic_factor(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    length(model) == 1 && return one(T)
    return PT_property(model,p,T,z,phase,threaded,vol0,VT_thermodynamic_factor)
end

#first derivative order properties
export entropy, internal_energy, enthalpy, gibbs_free_energy, helmholtz_free_energy
export entropy_res, internal_energy_res, enthalpy_res, gibbs_free_energy_res, helmholtz_free_energy_res
export gibbs_energy,helmholtz_energy,gibbs_energy_res,helmholtz_energy_res
export mass_enthalpy,mass_entropy,mass_internal_energy,mass_gibbs_energy,mass_gibbs_free_energy,mass_helmholtz_energy,mass_helmholtz_free_energy
#second derivative order properties
export isochoric_heat_capacity, isobaric_heat_capacity,adiabatic_index
export mass_isobaric_heat_capacity,mass_isochoric_heat_capacity
export isothermal_compressibility, isentropic_compressibility, speed_of_sound
export isobaric_expansivity, joule_thomson_coefficient, inversion_temperature
#higher derivative order properties
export fundamental_derivative_of_gas_dynamics
#volume properties
export mass_density,molar_density, compressibility_factor
#molar gradient properties
export chemical_potential, activity_coefficient, activity, aqueous_activity, fugacity_coefficient,reference_chemical_potential,reference_chemical_potential_type
export chemical_potential_res
export mixing, excess, gibbs_solvation, partial_property
export identify_phase
export thermodynamic_factor

"""
    PT0

Module that stores Clapeyron properties in (total) volume-temperature basis.

All bulk properties have the following form:

```julia
property(model,p,T,z,phase=:unknown, threaded=true, vol0=nothing)
```

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

only a volume solver is done to get the volume from the corresponding pressure-temperature pair.
For a module that does a flash to check if there are more than one phase use the `PT` module instead.
"""
module PT0
    import Clapeyron
    for prop in Clapeyron.CLAPEYRON_PROPS
        @eval begin
            function $prop(model, p, T, z = Clapeyron.SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
                return Clapeyron.$prop(model,p,T,z;phase,threaded,vol0)
            end
        end
    end

    function flash(model,p,T,z = Clapeyron.SA[1.0],args...;kwargs...)
        return Clapeyron.tp_flash2(model,p,T,z,args...;kwargs...)
    end
end #module


function PT_property_withflash(model,p,T,z,phase,f::F) where {F}
    if f == pressure
        return p
    elseif f == temperature
        return T
    end

    if z isa Number
        return PT_property_withflash(model,p,T,SA[z],phase,f)
    end
    if !is_unknown(phase)
        V = volume(model, p, T, z; phase)
        return PT_property(model,p,T,z,phase,volume(model, p, T, z; phase),f) 
    else
        res = tp_flash2(model,p,T,z)
        ff = PT_to_VT(f)
        return ff(model,res)
    end
end

"""
    PT

Module that stores Clapeyron properties in pressure-temperature basis.

All bulk properties have the following form:

```julia
property(model,p,T,z,phase=:unknown)
```

If no `phase` argument is passed as an input, a pressure-temperature flash is done to check if the input pair corresponds to one or more phases.
To evaluate the property directly in the P-T base, use the `PT0` module instead.
"""
module PT
    import Clapeyron
    for prop in Clapeyron.CLAPEYRON_PROPS
        @eval begin
            function $prop(model, p, T, z = Clapeyron.SA[1.]; phase=:unknown)
                return Clapeyron.PT_property_withflash(model,p,T,z,phase,Clapeyron.$prop)
            end
        end
    end

    function flash(model,p,T,z = Clapeyron.SA[1.0],args...;kwargs...)
        return Clapeyron.tp_flash2(model,p,T,z,args...;kwargs...)
    end
end #module

"""
    supports_lever_rule(::f)::Bool

Returns `true` if the input property function can be used to describe multiphase mixtures using the lever rule:
```
f(a)/f(b) = (f - f(b))/f(a) - f(b))
```
"""
supports_lever_rule(f) = false

for prop in [:volume, :pressure, :entropy, :internal_energy, :enthalpy, :gibbs_free_energy, :helmholtz_free_energy,
    :entropy_res, :internal_energy_res, :enthalpy_res, :gibbs_free_energy_res, :helmholtz_free_energy_res,
   #volume :properties
    :mass_density,:molar_density,
    :mass_enthalpy,:mass_entropy,:mass_internal_energy,:mass_gibbs_free_energy,:mass_helmholtz_free_energy]
    @eval begin
        supports_lever_rule(::typeof($prop)) = true
    end
end

function spec_to_vt end

for prop in CLAPEYRON_PROPS
    prop in CLAPEYRON_PROP_ALIASES && continue
    VT_prop = VT_symbol(prop)
    @eval begin
        function spec_to_vt(model,V,T,z,spec::typeof($prop))
            VT0.$prop(model,V,T,z)
        end

        PT_to_VT(x::typeof($prop)) = $VT_prop
    end
end