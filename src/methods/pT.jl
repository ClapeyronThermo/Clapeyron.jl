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
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_entropy(model, V, T, z)
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
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_entropy_res(model, V, T, z)
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
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_chemical_potential(model,V,T,z)
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
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_chemical_potential_res(model,V,T,z)
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
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_internal_energy(model,V,T,z)
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
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_internal_energy_res(model,V,T,z)
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
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_enthalpy(model,V,T,z)
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
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_enthalpy_res(model,V,T,z)
end

"""
    gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the gibbs free energy, defined as:

```julia
G = A + p*V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_gibbs_free_energy(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    A = eos(model,V,T,z)
    return A + V*p
    #return A - V*∂A∂V
    #return VT_gibbs_free_energy(model,V,T,z)
end

"""
    gibbs_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the residual gibbs free energy, defined as:

```julia
G = Ar - V*∂Ar/∂V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_gibbs_free_energy_res(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function gibbs_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_gibbs_free_energy_res(model,V,T,z)
end

"""
    helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the helmholtz free energy, defined as:

```julia
A = eos(model,V(p),T,z)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_helmholtz_free_energy(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_helmholtz_free_energy(model,V,T,z)
end

"""
    helmholtz_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J]`

Calculates the residual helmholtz free energy, defined as:

```julia
A = eos_res(model,V(p),T,z)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_helmholtz_free_energy_res(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function helmholtz_free_energy_res(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_helmholtz_free_energy_res(model,V,T,z)
end

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
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).
"""
function isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_isochoric_heat_capacity(model,V,T,z)
end

"""
    isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Default units: `[J/K]`

Calculates the isobaric heat capacity, defined as:

```julia
Cp =  -T*(∂²A/∂T² - (∂²A/∂V∂T)^2 / ∂²A/∂V²)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isobaric_heat_capacity(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).

"""
function isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_isobaric_heat_capacity(model,V,T,z)
end
"""
    isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[Pa^-1]`

Calculates the isothermal compressibility, defined as:

```julia
κT =  (V*∂p/∂V)^-1
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isothermal_compressibility(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_isothermal_compressibility(model,V,T,z)
end

"""
    isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[Pa^-1]`

Calculates the isentropic compressibility, defined as:

```julia
κS =  (V*( ∂²A/∂V² - ∂²A/∂V∂T^2 / ∂²A/∂T² ))^-1
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isentropic_compressibility(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).

"""
function isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_isentropic_compressibility(model,V,T,z)
end

"""
    speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[m/s]`

Calculates the speed of sound, defined as:

```julia
c =  V * √(∂²A/∂V² - ∂²A/∂V∂T^2 / ∂²A/∂T²)/Mr)
```
Where `Mr` is the molecular weight of the model at the input composition.

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_speed_of_sound(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).

"""
function speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_speed_of_sound(model,V,T,z)
end

"""
    isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[K^-1]`

Calculates the isobaric expansivity, defined as:

```julia
α =  -∂²A/∂V∂T / (V*∂²A/∂V²)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isobaric_expansivity(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_isobaric_expansivity(model,V,T,z)
end

"""
    joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

default units: `[K/Pa]`

Calculates the joule thomson coefficient, defined as:

```julia
μⱼₜ =  -(∂²A/∂V∂T - ∂²A/∂V² * ((T*∂²A/∂T² + V*∂²A/∂V∂T) / (T*∂²A/∂V∂T + V*∂²A/∂V²)))^-1
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_joule_thomson_coefficient(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).

"""
function joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_joule_thomson_coefficient(model,V,T,z)
end

"""
    fugacity_coefficient(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)

Calculates the fugacity coefficient φᵢ, defined as:

```julia
log(φᵢ) =  μresᵢ/RT - log(Z)
```
Where `μresᵢ` is the vector of residual chemical potentials and `Z` is the compressibility factor.

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_fugacity_coefficient(model,V,T,z)`.

The keywords `phase`, `threaded` and `vol0` are passed to the [`Clapeyron.volume`](@ref) solver.
"""
function fugacity_coefficient(model::EoSModel,p,T,z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_fugacity_coefficient(model,V,T,z)
end


function fugacity_coefficient!(φ,model::EoSModel,p,T,z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    V = volume(model, p, T, z; phase, threaded, vol0)
    φ = VT_chemical_potential_res!(φ,model,V,T,z)
    R̄ = Rgas(model)
    Z = p*V/R̄/T/sum(z)
    φ ./= (R̄*T)
    φ .= exp.(φ)
    φ ./= Z
    return φ
end

function activity_coefficient(model::EoSModel,p,T,z=SA[1.]; phase=:unknown, threaded=true, vol0=nothing)
    pure   = split_model(model)
    μ_mixt = chemical_potential(model, p, T, z; phase, threaded, vol0)
    μ_pure = gibbs_free_energy.(pure, p, T; phase, threaded, vol0)
    R̄ = Rgas(model)
    return exp.((μ_mixt .- μ_pure) ./ R̄ ./ T) ./z
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
    V = volume(model, p, T, z; phase, threaded, vol0)
    R̄ = Rgas(model)
    return p*V/(sum(z)*R̄*T)
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
ρₙ =  ∑nᵢ/V
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
ρₙ =  Mr/V
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

function partial_property(model::EoSModel, p, T, z, property::ℜ; phase=:unknown, threaded=true, vol0=nothing) where {ℜ}
    V = volume(model, p, T, z; phase, threaded, vol0)
    return VT_partial_property(model,V,T,z,property)
end

#special dispatch for volume here
function VT_partial_property(model::EoSModel, V, T, z, ::typeof(volume))
    _,dpdv = p∂p∂V(model,V,T,z)
    dpdni = VT_partial_property(model, V, T, z, pressure)
    return -dpdni ./ dpdv
end

#first derivative order properties
export entropy, internal_energy, enthalpy, gibbs_free_energy, helmholtz_free_energy
export entropy_res, internal_energy_res, enthalpy_res, gibbs_free_energy_res, helmholtz_free_energy_res
#second derivative order properties
export isochoric_heat_capacity, isobaric_heat_capacity
export isothermal_compressibility, isentropic_compressibility, speed_of_sound
export isobaric_expansivity, joule_thomson_coefficient, inversion_temperature
#volume properties
export mass_density,molar_density, compressibility_factor
#molar gradient properties
export chemical_potential, activity_coefficient, fugacity_coefficient
export mixing, excess, gibbs_solvation