function PT_property(model,p,T,z,phase,threaded,vol0,f::F,::Val{UseP}) where {F,UseP}
    if z isa Number
        return PT_property(model,p,T,SA[z],phase,threaded,vol0,f,Val{UseP}())
    end
    V = volume(model, p, T, z; phase, threaded, vol0)
    if UseP
        return f(model,V,T,z,p)
    else
        return f(model,V,T,z)
    end
end

function PT_property(model,p,T,z,phase,threaded,vol0,f::F) where F
    PT_property(model,p,T,z,phase,threaded,vol0,f,Val{false}())
end

"""
    entropy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J/K]`

Calculates entropy, defined as:

```julia
S = -∂A/∂T
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_entropy(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function entropy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    return PT_property(model,p,T,z,phase,threaded,vol0,VT_entropy)
end

"""
    entropy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J/K]`

Calculates residual entropy, defined as:

```julia
S = -∂Ares/∂T
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_entropy_res(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function entropy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    return PT_property(model,p,T,z,phase,threaded,vol0,VT_entropy_res)
end

"""
    chemical_potential(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J/mol]`

Calculates the chemical potential, defined as:

```julia
μᵢ = ∂A/∂nᵢ
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_chemical_potential(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function chemical_potential(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_chemical_potential)
end

"""
    chemical_potential_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J/mol]`

Calculates the residual chemical potential, defined as:

```julia
μresᵢ = ∂Ares/∂nᵢ
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_chemical_potential_res(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function chemical_potential_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_chemical_potential_res)
end

"""
    internal_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the internal energy, defined as:

```julia
U = A - T * ∂A/∂T
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_internal_energy(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function internal_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_internal_energy)
end

"""
    internal_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the residual internal energy, defined as:

```julia
U = Ar - T * ∂Ar/∂T
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_internal_energy_res(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function internal_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_internal_energy_res_res)
end

"""
    enthalpy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the enthalpy, defined as:

```julia
H = A - T * ∂A/∂T - V * ∂A/∂V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_enthalpy(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function enthalpy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_enthalpy)
end

"""
    enthalpy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the residual enthalpy, defined as:

```julia
H = Ar - T * ∂Ar/∂T - V * ∂Ar/∂V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_enthalpy_res(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function enthalpy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_enthalpy_res)
end

"""
    gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    gibbs_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the gibbs free energy, defined as:

```julia
G = A + p*V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_gibbs_free_energy(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_gibbs_free_energy,Val{true}())
end

"""
    gibbs_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    gibbs_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the residual gibbs free energy, defined as:

```julia
G = Ar - V*∂Ar/∂V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_gibbs_free_energy_res(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function gibbs_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_gibbs_free_energy_res)
end

"""
    helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    helmholtz_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the helmholtz free energy, defined as:

```julia
A = eos(model,V(p),T,z)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_helmholtz_free_energy(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_helmholtz_free_energy)
end

"""
    helmholtz_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    helmholtz_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the residual helmholtz free energy, defined as:

```julia
A = eos_res(model,V(p),T,z)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_helmholtz_free_energy_res(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function helmholtz_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_helmholtz_free_energy_res)
end

const helmholtz_energy = helmholtz_free_energy
const helmholtz_energy_res = helmholtz_free_energy_res 
const gibbs_energy = gibbs_free_energy
const gibbs_energy_res = gibbs_free_energy_res

"""
    isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J/K]`

Calculates the isochoric heat capacity, defined as:

```julia
Cv = -T * ∂²A/∂T²
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isochoric_heat_capacity(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel = ReidIdeal)`).
"""
function isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isochoric_heat_capacity)
end

"""
    isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J/K]`

Calculates the isobaric heat capacity, defined as:

```julia
Cp = -T*(∂²A/∂T² - (∂²A/∂V∂T)^2 / ∂²A/∂V²)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isobaric_heat_capacity(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel = ReidIdeal)`).

"""
function isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isobaric_heat_capacity)
end

"""
    adiabatic_index(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J/K]`

Calculates the isobaric heat capacity, defined as:

```julia
γ = Cp/Cv
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_adiabatic_index(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel = ReidIdeal)`).

"""
function adiabatic_index(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_adiabatic_index)
end

"""
    isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[Pa^-1]`

Calculates the isothermal compressibility, defined as:

```julia
κT = -(V*∂p/∂V)^-1
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isothermal_compressibility(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isothermal_compressibility)
end

"""
    isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[Pa^-1]`

Calculates the isentropic compressibility, defined as:

```julia
κS = (V*( ∂²A/∂V² - ∂²A/∂V∂T^2 / ∂²A/∂T² ))^-1
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isentropic_compressibility(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel = ReidIdeal)`).

"""
function isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isentropic_compressibility)
end

"""
    speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[m/s]`

Calculates the speed of sound, defined as:

```julia
c = V * √(∂²A/∂V² - ∂²A/∂V∂T^2 / ∂²A/∂T²)/Mr)
```
Where `Mr` is the molecular weight of the model at the input composition.

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_speed_of_sound(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel = ReidIdeal)`).

"""
function speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_speed_of_sound)
end

"""
    isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[K^-1]`

Calculates the isobaric expansivity, defined as:

```julia
α = -∂²A/∂V∂T / (V*∂²A/∂V²)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isobaric_expansivity(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_isobaric_expansivity)
end

"""
    joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[K/Pa]`

Calculates the joule thomson coefficient, defined as:

```julia
μⱼₜ = -(∂²A/∂V∂T - ∂²A/∂V² * ((T*∂²A/∂T² + V*∂²A/∂V∂T) / (T*∂²A/∂V∂T + V*∂²A/∂V²)))^-1
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_joule_thomson_coefficient(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel = ReidIdeal)`).

"""
function joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_joule_thomson_coefficient)
end

"""
    identify_phase(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)::Symbol


Returns the phase of a fluid at the conditions specified by `V`, `T` and `z`.
Uses the phase identification parameter criteria from `Clapeyron.pip`

returns `:liquid` if the phase is liquid (or liquid-like), `:vapour` if the phase is vapour (or vapour-like), and `:unknown` if the calculation of the phase identification parameter failed.

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_enthalpy(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function identify_phase(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    #TODO: what do we do with composite models here?
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_identify_phase(model,V,T,z)
end

"""
    fundamental_derivative_of_gas_dynamics(model::EoSModel, p, T, z=SA[1.]; phase=:gas, threaded=true, vol0=nothing)::Symbol

Calculates the fundamental derivative of gas dynamics.

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_enthalpy(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function fundamental_derivative_of_gas_dynamics(model::EoSModel, p, T, z=SA[1.]; phase=:gas, threaded=true, vol0=nothing)
    PT_property(model,p,T,z,phase,threaded,vol0,VT_fundamental_derivative_of_gas_dynamics)
end

"""
    fugacity_coefficient(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Calculates the fugacity coefficient φᵢ, defined as:

```julia
log(φᵢ) = μresᵢ/RT - log(Z)
```
Where `μresᵢ` is the vector of residual chemical potentials and `Z` is the compressibility factor.

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_fugacity_coefficient(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function fugacity_coefficient(model::EoSModel,p,T,z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
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
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_fugacity_coefficient(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

If the `μ_ref` keyword argument is not provided, the `reference` keyword is used to specify the reference chemical potential..
"""
function activity_coefficient(model::EoSModel,p,T,z=SA[1.];
                            μ_ref = nothing,
                            reference = :pure,
                            phase=:unknown,
                            threaded=true,
                            vol0=nothing)
    if model isa ActivityModel
        return activity_coefficient(model,p,T,z)
    end
    if μ_ref == nothing
        return activity_coefficient_impl(model,p,T,z,reference_chemical_potential(model,p,T,reference;phase,threaded),reference,phase,threaded,vol0)
    else
        return activity_coefficient_impl(model,p,T,z,μ_ref,reference,phase,threaded,vol0)
    end
end

function activity_coefficient_impl(model,p,T,z,μ_ref,reference,phase,threaded,vol0)
    R̄ = Rgas(model)
    μ_mixt = chemical_potential(model, p, T, z; phase, threaded, vol0)
    return sum(z) .* exp.((μ_mixt .- μ_ref) ./ R̄ ./ T) ./z
end

"""
    activity(model::EoSModel,p,T,z=SA[1.0];reference = :pure, phase=:unknown, threaded=true, vol0=nothing)

Calculates the activity, defined as:
```julia
log(a) = (μ_mixt - μ_ref) / R̄ / T
```
where `μ_mixt` is the chemical potential of the mixture and `μ_ref` is the reference chemical potential for the model at `p`,`T` conditions, calculated via [`Clapeyron.reference_chemical_potential`](@ref).
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_fugacity_coefficient(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

If the `μ_ref` keyword argument is not provided, the `reference` keyword is used to specify the reference chemical potential..
"""
function activity(model::EoSModel,p,T,z=SA[1.];
                μ_ref = nothing,
                reference = :pure,
                phase=:unknown,
                threaded=true,
                vol0=nothing)
    if model isa ActivityModel
        return activity(model,p,T,z)
    end
    if μ_ref == nothing
        return activity_impl(model,p,T,z,reference_chemical_potential(model,p,T,reference;phase,threaded),reference,phase,threaded,vol0)
    else
        return activity_impl(model,p,T,z,μ_ref,reference,phase,threaded,vol0)
    end
end

function activity_impl(model,p,T,z,μ_ref,reference,phase,threaded,vol0)
    R̄ = Rgas(model)
    μ_mixt = chemical_potential(model, p, T, z; phase, threaded, vol0)
    return exp.((μ_mixt .- μ_ref) ./ R̄ ./ T)
end

function find_hydronium_index(model)
    idx = findfirst(isequal("hydronium"),model.components)
    idx == nothing && return 0
    return idx
end

function find_hydroxide_index(model)
    idx = findfirst(isequal("hydroxide"),model.components)
    idx == nothing && return 0
    return idx
end

function find_water_indx(model)
    idx = findfirst(isequal("water"),model.components)
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

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_fugacity_coefficient(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
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

Returns a reference chemical potential. used in calculation of `activity` and actitivy_coefficient. there are two available references:
- `:pure`: the reference potential is a pure component at specified `T`, `p` and `phase`
- `:aqueous`: the chemical potential of the pure components at specified `T`, `p` and `phase`
- `:sat_pure_T`:  the reference potential is the pure saturated liquid phase at specified `T`.
- `:zero`: the reference potential is equal to zero for all components (used for `ActivityModel`)
The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function reference_chemical_potential(model::EoSModel,p,T,reference = reference_chemical_potential_type(model); phase=:unknown, threaded=true, vol0=nothing)
    if reference == :pure
        pure = split_model.(model)
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
        pure = split_model.(model)
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
    compressibility_factor(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Calculates the compressibility factor `Z`, defined as:

```julia
Z = p*V(p)/R*T
```
The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function compressibility_factor(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    #this property only depends on the implementation of volume_impl.
    V = volume(model,p,T,z;phase,threaded,vol0)
    return p*V/(sum(z)*Rgas(model)*T)
end

function inversion_temperature(model::EoSModel, p, z=SA[1.0]; phase=:unknown, threaded=true, vol0=nothing)
    T0 = 6.75*T_scale(model,z)
    μⱼₜ(T) = joule_thomson_coefficient(model, p, T, z; phase, threaded, vol0)
    return Roots.find_zero(μⱼₜ,T0)
end

"""
    molar_density(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[mol/m^3]`

Calculates the molar density, defined as:

```julia
ρₙ = ∑nᵢ/V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_molar_density(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function molar_density(model::EoSModel,p,T,z=SA[1.0];phase=:unknown, threaded=true, vol0=nothing)
     V = volume(model, p, T, z; phase, threaded, vol0)
     return VT_molar_density(model,V,T,z)
end

"""
    mass_density(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true)

default units: `[kg/m^3]`

Calculates the mass density, defined as:

```julia
ρₙ = Mr/V
```
Where `Mr` is the molecular weight of the model at the input composition.

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_mass_density(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function mass_density(model::EoSModel, p, T, z=SA[1.0]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_mass_density(model,V,T,z)
end

"""
    mixing(model::EoSModel, p, T, z=SA[1.], property; phase=:unknown, threaded=true, vol0=nothing)

Calculates the mixing function for a specified property as:

```julia
f_mix = f(p,T,z) - ∑zᵢ*f_pureᵢ(p,T)
```
The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function mixing(model::EoSModel, p, T, z, property::ℜ; phase=:unknown, threaded=true, vol0=nothing) where {ℜ}
    pure = split_model(model)
    TT = typeof(p+T+first(z))
    mix_prop  = property(model, p, T, z; phase, threaded, vol0)
    for i in 1:length(z)
        mix_prop -= z[i]*property(pure[i], p, T; phase, threaded, vol0)
    end
    return mix_prop::TT
end

function excess(model::EoSModel, p, T, z, property; phase=:unknown, threaded=true, vol0=nothing)
    mixing(model, p, T, z, property; phase, threaded, vol0)
end

function excess(model::EoSModel, p, T, z, ::typeof(entropy); phase=:unknown, threaded=true, vol0=nothing)
    TT = typeof(p+T+first(z))
    pure = split_model(model)
    s_mix = entropy_res(model, p, T, z; phase, threaded, vol0)
    for i in 1:length(z)
        s_mix -= z[i]*entropy_res(pure[i], p, T; phase, threaded, vol0)
    end
    #s_pure = entropy_res.(pure,p,T)
    return s_mix::TT
end

function excess(model::EoSModel, p, T, z, ::typeof(gibbs_free_energy); phase=:unknown, threaded=true, vol0=nothing)
    TT = typeof(p+T+first(z))
    pure = split_model(model)
    g_mix = gibbs_free_energy(model, p, T, z; phase, threaded, vol0)
    log∑z = log(sum(z))
    R̄ = Rgas(model)
    for i in 1:length(z)
        lnxi = R̄*T*(log(z[i]) - log∑z)
        g_mix -= z[i]*(gibbs_free_energy(pure[i], p, T; phase, threaded, vol0) + lnxi)
    end

    return g_mix::TT
end


"""
    gibbs_solvation(model::EoSModel, T; threaded=true, vol0=(nothing,nothing))

Calculates the solvation free energy as:

```julia
g_solv = -R̄*T*log(K)
```
where the first component is the solvent and second is the solute.
"""
function gibbs_solvation(model::EoSModel, T; threaded=true, vol0=(nothing,nothing))
    binary_component_check(gibbs_solvation, model)
    pure = split_model(model)
    z = [1.0,1e-30]

    p,v_l,v_v = saturation_pressure(pure[1],T)

    φ_l = fugacity_coefficient(model, p, T, z; phase=:l, threaded, vol0=vol0[1])
    φ_v = fugacity_coefficient(model, p, T, z; phase=:v, threaded, vol0=vol0[2])

    K = φ_v[2]*v_v/φ_l[2]/v_l
    R̄ = Rgas(model)
    return -R̄*T*log(K)
end

"""
    partial_property(model::EoSModel, p, T, z, property::X; phase=:unknown, threaded=true, vol0=nothing) where {X}

Calculate the partial molar property of a mixture at specified temperature, pressure, mol amounts, and extensive property of interest.
The equality `sum(z .* partial_property(model,p,T,z,property) - property(model,p,T,z))` should hold.
    
The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function partial_property(model::EoSModel, p, T, z, property::ℜ; phase=:unknown, threaded=true, vol0=nothing) where {ℜ}
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
    xv(∂V) = VT_prop(model,∂V,T,z)
    ∂x∂V = Solvers.derivative(xv,V)
    _,∂p∂V = p∂p∂V(model,V,T,z)
    return ∂x∂nᵢ .- ∂x∂V .* ∂p∂nᵢ ./ ∂p∂V
end

#default
PT_to_VT(x) = x

for (PTprop,VTprop) in [
    (:entropy,:VT_entropy),
    (:enthalpy,:VT_enthalpy),
    (:internal_energy,:VT_internal_energy),
    (:gibbs_free_energy,:VT_gibbs_free_energy),
    (:helmholtz_free_energy,:VT_helmholtz_free_energy),
    (:isochoric_heat_capacity,:VT_isochoric_heat_capacity),
    (:isobaric_heat_capacity,:VT_isobaric_heat_capacity)
    ]
    @eval begin
        PT_to_VT(::typeof($PTprop)) = $VTprop
    end
end

#first derivative order properties
export entropy, internal_energy, enthalpy, gibbs_free_energy, helmholtz_free_energy
export entropy_res, internal_energy_res, enthalpy_res, gibbs_free_energy_res, helmholtz_free_energy_res
export gibbs_energy,helmholtz_energy,gibbs_energy_res,helmholtz_energy_res

#second derivative order properties
export isochoric_heat_capacity, isobaric_heat_capacity,adiabatic_index
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

module PT
    #first derivative order properties
    using Clapeyron: entropy, internal_energy, enthalpy, gibbs_free_energy, helmholtz_free_energy
    using Clapeyron: entropy_res, internal_energy_res, enthalpy_res, gibbs_free_energy_res, helmholtz_free_energy_res
    #second derivative order properties
    using Clapeyron: isochoric_heat_capacity, isobaric_heat_capacity,adiabatic_index
    using Clapeyron: isothermal_compressibility, isentropic_compressibility, speed_of_sound
    using Clapeyron: isobaric_expansivity, joule_thomson_coefficient, inversion_temperature
    #higher derivative order properties
    using Clapeyron: fundamental_derivative_of_gas_dynamics
    #volume properties
    using Clapeyron: mass_density,molar_density, compressibility_factor
    using Clapeyron: identify_phase
    import Clapeyron
    pressure(model, p, T, z=Clapeyron.SA[1.]; phase=:unknown, threaded=true, vol0=nothing) = p
    temperature(model, p, T, z=Clapeyron.SA[1.]; phase=:unknown, threaded=true, vol0=nothing) = T
    function flash(model,p,T,z = Clapeyron.SA[1.0],args...;kwargs...)
        return Clapeyron.tp_flash2(model,p,T,z,args...;kwargs...)
    end
end

"""
    supports_lever_rule(::f)::Bool

returns `true` if the input property function can be used to describe multiphase mixtures using the lever rule:
```
f(a)/f(b) = (f - f(b))/f(a) - f(b))
```
"""
supports_lever_rule(f) = false

for prop in [:volume, :pressure, :entropy, :internal_energy, :enthalpy, :gibbs_free_energy, :helmholtz_free_energy,
    :entropy_res, :internal_energy_res, :enthalpy_res, :gibbs_free_energy_res, :helmholtz_free_energy_res,
   #volume :properties
    :mass_density,:molar_density]
    @eval begin
        supports_lever_rule(::typeof($prop)) = true
    end
end
