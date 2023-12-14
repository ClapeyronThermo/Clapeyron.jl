"""
    entropy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J/K]`

Calculates entropy, defined as:

```julia
S = -∂A/∂T
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_entropy(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function entropy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_entropy(model,V,T,z)
end

"""
    entropy_res(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J/K]`

Calculates residual entropy, defined as:

```julia
S = -∂Ares/∂T
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_entropy_res(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function entropy_res(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_entropy_res(model,V,T,z)
end

"""
    chemical_potential(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J/mol]`    

Calculates the chemical potential, defined as:

```julia
μᵢ = ∂A/∂nᵢ
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_chemical_potential(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function chemical_potential(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_chemical_potential(model,V,T,z)
end

"""
    chemical_potential_res(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J/mol]`

Calculates the residual chemical potential, defined as:

```julia
μresᵢ = ∂Ares/∂nᵢ
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_chemical_potential_res(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function chemical_potential_res(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_chemical_potential_res(model,V,T,z)
end

"""
    internal_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J]`

Calculates the internal energy, defined as:

```julia
U = A - T * ∂A/∂T
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_internal_energy(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function internal_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_internal_energy(model,V,T,z)
end

"""
    enthalpy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J]`
     
Calculates the enthalpy, defined as:

```julia
H = A - T * ∂A/∂T - V * ∂A/∂V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_enthalpy(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function enthalpy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_enthalpy(model,V,T,z)
end

"""
    gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J]`
     
Calculates the gibbs free energy, defined as:

```julia
G = A - V * ∂A/∂V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `VT_gibbs_free_energy(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function gibbs_free_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    A = eos(model,V,T,z)
    return A + V*p
    #return A - V*∂A∂V
    #return VT_gibbs_free_energy(model,V,T,z)
end

"""
    helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J]`
     
Calculates the helmholtz free energy, defined as:

```julia
A = eos(model,V(p),T,z)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and calculates the property via `eos(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function helmholtz_free_energy(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_helmholtz_free_energy(model,V,T,z)
end

"""
    isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J/K]`

Calculates the isochoric heat capacity, defined as:

```julia
Cv = -T * ∂²A/∂T²
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isochoric_heat_capacity(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).
"""
function isochoric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isochoric_heat_capacity(model,V,T,z)
end

"""
    isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Default units: `[J/K]`

Calculates the isobaric heat capacity, defined as:

```julia
Cp =  -T*(∂²A/∂T² - (∂²A/∂V∂T)^2 / ∂²A/∂V²)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and
calculates the property via `VT_isobaric_heat_capacity(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).

"""
function isobaric_heat_capacity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isobaric_heat_capacity(model,V,T,z)
end
"""
    isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

default units: `[Pa^-1]`

Calculates the isothermal compressibility, defined as:

```julia
κT =  (V*∂p/∂V)^-1
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and 
calculates the property via `VT_isothermal_compressibility(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function isothermal_compressibility(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isothermal_compressibility(model,V,T,z)
end

"""
    isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

default units: `[Pa^-1]`    

Calculates the isentropic compressibility, defined as:

```julia
κS =  (V*( ∂²A/∂V² - ∂²A/∂V∂T^2 / ∂²A/∂T² ))^-1
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and 
calculates the property via `VT_isentropic_compressibility(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).

"""
function isentropic_compressibility(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isentropic_compressibility(model,V,T,z)
end

"""
    speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

default units: `[m/s]`

Calculates the speed of sound, defined as:

```julia
c =  V * √(∂²A/∂V² - ∂²A/∂V∂T^2 / ∂²A/∂T²)/Mr)
```
Where `Mr` is the molecular weight of the model at the input composition.

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and 
calculates the property via `VT_speed_of_sound(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).

"""
function speed_of_sound(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_speed_of_sound(model,V,T,z)
end

"""
    isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

default units: `[K^-1]`

Calculates the isobaric expansivity, defined as:

```julia
α =  -∂²A/∂V∂T / (V*∂²A/∂V²)
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and 
calculates the property via `VT_isobaric_expansivity(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function isobaric_expansivity(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_isobaric_expansivity(model,V,T,z)
end

"""
    joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

default units: `[K/Pa]`

Calculates the joule thomson coefficient, defined as:

```julia
μⱼₜ =  -(∂²A/∂V∂T - ∂²A/∂V² * ((T*∂²A/∂T² + V*∂²A/∂V∂T) / (T*∂²A/∂V∂T + V*∂²A/∂V²)))^-1
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and 
calculates the property via `VT_joule_thomson_coefficient(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.

!!! warning "Accurate ideal model required"
    This property requires at least second order ideal model temperature derivatives. If you are computing these properties, consider using a different ideal model than the `BasicIdeal` default (e.g. `EoS(["species"];idealmodel=ReidIdeal)`).

"""
function joule_thomson_coefficient(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    return VT_joule_thomson_coefficient(model,V,T,z)
end

"""
    fugacity_coefficient(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)
     
Calculates the fugacity coefficient φᵢ, defined as:

```julia
log(φᵢ) =  μresᵢ/RT - log(Z)
```
Where `μresᵢ` is the vector of residual chemical potentials and `Z` is the compressibility factor.

The keywords `phase` and `threaded` are passed to [`Clapeyron.volume`](@ref).
"""
function fugacity_coefficient(model::EoSModel,p,T,z=SA[1.]; phase = :unknown, threaded=true)
    V = volume(model,p,T,z;phase,threaded)
    μ_res = VT_chemical_potential_res(model,V,T,z)
    φ = μ_res
    R̄ = Rgas(model)
    Z = p*V/R̄/T/sum(z)
    return exp.(μ_res ./ R̄ ./ T) ./ Z
end

function fugacity_coefficient!(φ,model::EoSModel,p,T,z=SA[1.]; phase = :unknown, threaded=true)
    V = volume(model,p,T,z;phase,threaded)
    φ = VT_chemical_potential_res!(φ,model,V,T,z)
    R̄ = Rgas(model)
    Z = p*V/R̄/T/sum(z)
    φ ./= (R̄*T)
    φ .= exp.(φ)
    φ ./= Z
    return φ
end

function activity_coefficient(model::EoSModel,p,T,z=SA[1.]; phase = :unknown, threaded=true)
    pure   = split_model(model)
    μ_mixt = chemical_potential(model,p,T,z;phase,threaded)
    μ_pure = gibbs_free_energy.(pure,p,T;phase,threaded)
    R̄ = Rgas(model)
    return exp.((μ_mixt .- μ_pure) ./ R̄ ./ T) ./z
end

"""
    compressibility_factor(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

Calculates the compressibility factor `Z`, defined as:

```julia
Z = p*V(p)/R*T
```
The keywords `phase` and `threaded` are passed to [`Clapeyron.volume`](@ref).
"""
function compressibility_factor(model::EoSModel, p, T, z=SA[1.]; phase = :unknown,threaded=true)
    V = volume(model, p, T, z; phase=phase, threaded=threaded)
    R̄ = Rgas(model)
    return p*V/(sum(z)*R̄*T)
end

function inversion_temperature(model::EoSModel, p, z=SA[1.0])
    T0 = 6.75*T_scale(model,z)
    μⱼₜ(T) = joule_thomson_coefficient(model,p,T,z)
    return Roots.find_zero(μⱼₜ,T0)
end

"""
    molar_density(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)

default units: `[mol/m^3]`

Calculates the molar density, defined as:

```julia
ρₙ =  ∑nᵢ/V
```
Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and 
calculates the property via `VT_molar_density(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function molar_density(model::EoSModel,p,T,z=SA[1.0];phase = :unknown,threaded=true)
     V = volume(model,p,T,z;phase=phase,threaded=threaded)
     return VT_molar_density(model,V,T,z)
end

"""
    mass_density(model::EoSModel, p, T, z=SA[1.]; phase = :unknown, threaded=true)
    
default units: `[kg/m^3]`

Calculates the molar density, defined as:

```julia
ρₙ =  Mr/V
```
Where `Mr` is the molecular weight of the model at the input composition.

Internally, it calls [`Clapeyron.volume`](@ref) to obtain `V` and 
calculates the property via `VT_mass_density(model,V,T,z)`.

The keywords `phase` and `threaded` are passed to the volume solver.
"""
function mass_density(model::EoSModel,p,T,z=SA[1.0];phase = :unknown,threaded=true)
    V = volume(model,p,T,z;phase=phase,threaded=threaded)
    return VT_mass_density(model,V,T,z)
end

"""
    mixing(model::EoSModel, p, T, z=SA[1.], property; phase = :unknown,threaded=true)

Calculates the mixing function for a specified property as:

```julia
f_mix = f(p,T,z) - ∑zᵢ*f_pureᵢ(p,T)
```
The keywords `phase` and `threaded` are passed to the volume solver.
"""
function mixing(model::EoSModel,p,T,z,property::ℜ;phase = :unknown,threaded=true) where {ℜ}
    pure = split_model(model)
    TT = typeof(p+T+first(z))
    mix_prop  = property(model,p,T,z;phase,threaded)
    for i in 1:length(z)
        mix_prop -= z[i]*property(pure[i],p,T;phase,threaded)
    end
    return mix_prop::TT
end

excess(model::EoSModel,p,T,z,property) = mixing(model::EoSModel,p,T,z,property)

function excess(model::EoSModel,p,T,z,::typeof(entropy))
    TT = typeof(p+T+first(z))
    pure = split_model(model)
    s_mix = entropy_res(model,p,T,z)
    for i in 1:length(z)
        s_mix -= z[i]*entropy_res(pure[i],p,T)
    end
    #s_pure = entropy_res.(pure,p,T)
    return s_mix::TT
end

function excess(model::EoSModel,p,T,z,::typeof(gibbs_free_energy))
    TT = typeof(p+T+first(z))
    pure = split_model(model)
    g_mix = gibbs_free_energy(model,p,T,z)
    log∑z = log(sum(z))
    v = volume(model,p,T,z)
    R̄ = Rgas(model)
    for i in 1:length(z)
        lnxi = R̄*T*(log(z[i]) - log∑z)
        g_mix -= z[i]*(gibbs_free_energy(pure[i],p,T) + lnxi)
    end

    return g_mix::TT
end


"""
    gibbs_solvation(model::EoSModel, T)

Calculates the solvation free energy as:

```julia
g_solv = -R̄*T*log(K)
```
where the first component is the solvent and second is the solute.
"""
function gibbs_solvation(model::EoSModel,T)
    binary_component_check(gibbs_solvation,model)
    pure = split_model(model)
    z = [1.0,1e-30]
    
    p,v_l,v_v = saturation_pressure(pure[1],T)
    
    φ_l = fugacity_coefficient(model,p,T,z;phase=:l)
    φ_v = fugacity_coefficient(model,p,T,z;phase=:v)
    
    K = φ_v[2]*v_v/φ_l[2]/v_l
    R̄ = Rgas(model)
    return -R̄*T*log(K)
end    

function partial_property(model::EoSModel,p,T,z,property::ℜ;phase = :unknown,threaded=true) where {ℜ}
    V = volume(model,p,T,z;phase,threaded)
    return VT_partial_property(model,V,T,z,property)
end

#special dispatch for volume here
function VT_partial_property(model::EoSModel,V,T,z,property::typeof(volume))
    _,dpdv = p∂p∂V(model,V,T,z)
    dpdni = VT_partial_property(model,V,T,z,pressure)
    return -dpdni ./ dpdv
end

export entropy, chemical_potential, internal_energy, enthalpy, gibbs_free_energy
export helmholtz_free_energy, isochoric_heat_capacity, isobaric_heat_capacity
export isothermal_compressibility, isentropic_compressibility, speed_of_sound
export isobaric_expansivity, joule_thomson_coefficient, compressibility_factor, inversion_temperature
export mass_density,molar_density, activity_coefficient, fugacity_coefficient, entropy_res
export mixing, excess, gibbs_solvation