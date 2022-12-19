"""
    pressure(model::EoSModel, V, T, z=SA[1.])

default units: `[Pa]`

Returns the pressure of the model at a given volume, temperature and composition, defined as:

```julia
p =  -∂A/∂V
```

"""
function pressure(model::EoSModel, V, T, z=SA[1.])
    return -∂f∂V(model,V,T,z)
end

function VT_entropy(model::EoSModel, V, T, z=SA[1.])
    return -∂f∂T(model,V,T,z)
end

function VT_entropy_res(model::EoSModel, V, T, z=SA[1.])
    fun(x) = eos_res(model,V,x,z)
    return -Solvers.derivative(fun,T)
end

function VT_internal_energy(model::EoSModel, V, T, z=SA[1.])
    dA, A = ∂f(model,V,T,z)
    ∂A∂V, ∂A∂T = dA
    return A - T*∂A∂T
end

function VT_enthalpy(model::EoSModel, V, T, z=SA[1.])
    dA, A = ∂f(model,V,T,z)
    ∂A∂V, ∂A∂T = dA
    return A - V*∂A∂V - T*∂A∂T
end

function VT_gibbs_free_energy(model::EoSModel, V, T, z=SA[1.])
    dA, A = ∂f(model,V,T,z)
    ∂A∂V, ∂A∂T = dA
    return A - V*∂A∂V
end

function VT_gibbs_free_energy_res(model::EoSModel, V, T, z=SA[1.])
    fun(x) = eos_res(model,x,T,z)
    Ar,∂A∂Vr = Solvers.f∂f(fun,V)
    return Ar - V*∂A∂Vr
end

function VT_helmholtz_free_energy(model::EoSModel, V, T, z=SA[1.])
    return eos(model,V,T,z)
end

function VT_helmholtz_free_energy_res(model::EoSModel, V, T, z=SA[1.])
    return eos_res(model,V,T,z)
end

function VT_isochoric_heat_capacity(model::EoSModel, V, T, z=SA[1.])
    A(x) = eos(model,V,x,z)
    ∂A∂T(x) = Solvers.derivative(A,x)
    ∂²A∂T²(x) = Solvers.derivative(∂A∂T,x)
    return -T*∂²A∂T²(T)
end

function VT_isobaric_heat_capacity(model::EoSModel, V, T, z=SA[1.])
    d²A = f_hess(model,V,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return -T*(∂²A∂T² - ∂²A∂V∂T^2/∂²A∂V²)
end

function VT_isothermal_compressibility(model::EoSModel, V, T, z=SA[1.])
    p0,∂p∂V = p∂p∂V(model,V,T,z)
    return -1/V/∂p∂V
end

function VT_isentropic_compressibility(model::EoSModel, V, T, z=SA[1.])
    d²A = f_hess(model,V,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return 1/V/(∂²A∂V²-∂²A∂V∂T^2/∂²A∂T²)
end

function VT_speed_of_sound(model::EoSModel, V, T, z=SA[1.])
    Mr = molecular_weight(model,z)
    d²A = f_hess(model,V,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return V*sqrt((∂²A∂V²-∂²A∂V∂T^2/∂²A∂T²)/Mr)
end

function VT_isobaric_expansivity(model::EoSModel, V, T, z=SA[1.])
    d²A = f_hess(model,V,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return -∂²A∂V∂T/(V*∂²A∂V²)
end

function VT_joule_thomson_coefficient(model::EoSModel, V, T, z=SA[1.])
    d²A = f_hess(model,V,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    return -(∂²A∂V∂T - ∂²A∂V²*((T*∂²A∂T² + V*∂²A∂V∂T) / (T*∂²A∂V∂T + V*∂²A∂V²)))^-1
end

"""
    second_virial_coefficient(model::EoSModel, T, z=SA[1.])

Default units: `[m^3]`

Calculates the second virial coefficient `B`, defined as:

```julia
B = lim(V->∞)[ V^2/RT *  (∂Aᵣ∂V + V*∂²Aᵣ∂V²) ]
```
where `Aᵣ` is the residual helmholtz energy.
"""
function second_virial_coefficient(model::EoSModel, T, z=SA[1.])
   return second_virial_coefficient_impl(model,T,z)
end

function second_virial_coefficient_impl(model::EoSModel,T , z = SA[1.0])
    TT = one(promote_type(eltype(z),typeof(1.0*T)))
    V = 1/sqrt(eps(TT))
    fAᵣ(x) = eos_res(model,x,T,z)
    Aᵣ,∂Aᵣ∂V,∂²Aᵣ∂V² = Solvers.f∂f∂2f(fAᵣ,V)
    return V^2/(R̄*T)*(∂Aᵣ∂V+V*∂²Aᵣ∂V²) #V*V/J * (J/V )
end


"""
    cross_second_virial(model,T,z)

Default units: `[m^3]`

Calculates the second cross virial coefficient (B₁₂) of a binary mixture, using the definition:

```julia
B̄ = x₁^2*B₁₁ + 2x₁x₂B₁₂ + x₂^2*B₂₂
B₁₂ = (B̄ - x₁^2*B₁₁ - x₂^2*B₂₂)/2x₁x₂
```


!!! info composition-dependent
    The second cross virial coefficient calculated from a equation of state can present a dependency on composition [1], but normally, experiments for obtaining the second virial coefficient are made by mixing the same volume of two gases. you can calculate B12 in this way by using (Clapeyron.equivol_cross_second_virial)[@ref]

## References
1. Jäger, A., Breitkopf, C., & Richter, M. (2021). The representation of cross second virial coefficients by multifluid mixture models and other equations of state. Industrial & Engineering Chemistry Research, 60(25), 9286–9295. [doi:10.1021/acs.iecr.1c01186](https://doi.org/10.1021/acs.iecr.1c01186)
"""
function cross_second_virial(model,T,z)
    B = second_virial_coefficient
    n = length(model)
    ∑z = sum(z)
    if n == 1
        return zero(T + first(z))
    elseif n == 2
        model1,model2 = split_model(model)
        B̄ = B(model,T,z)/∑z #1 mol
        B1,B2 = B(model1,T),B(model2,T) #1 mol by default
        #@show B1,B2,B̄
        x = z/∑z
        #B̄ = (B1*x1^2 + B2*x2^2 + B12*x1*x2)
        B12 = (B̄ - x[1]*x[1]*B1 - x[2]*x[2]*B2)/(2*x[1]*x[2])
        return B12*∑z
    else
        throw(error("cross_second_virial is only for models with 2 components. got a model with $n conponents"))
    end
end

"""
    equivol_cross_second_virial(model::EoSModel,T,p_exp = 200000.0)

calculates the second cross virial coefficient, by simulating the mixing of equal volumes of pure gas, at T,P conditions.
The equal volume of each pure gas sets an specific molar amount for each component. Details of the experiment can be found at [1].

## Example
```
model = SAFTVRQMie(["helium","neon"])
B12 = equivol_cross_second_virial(model,)

```
## References
1. Brewer, J., & Vaughn, G. W. (1969). Measurement and correlation of some interaction second virial coefficients from − 125° to 50°C. I. The Journal of Chemical Physics, 50(7), 2960–2968. [doi:10.1063/1.1671491](https://doi.org/10.1063/1.1671491)
"""
function equivol_cross_second_virial(model,T,p_exp = 200000.0)
    @assert length(model) == 2 "this function only works with binary models"
    #they do experiments at constant volume and temperature, so we are gonna need to calculate the mole fractions for that
    m1,m2 = split_model(model)
    B = second_virial_coefficient
    B11,B22 = B(m1,T),B(m2,T)
    v01,v02 = volume_virial(B11,p_exp,T),volume_virial(B22,p_exp,T)
    v1,v2 = volume(m1,p_exp,T,vol0 = v01),volume(m1,p_exp,T,vol0 = v02) #we do this in case there is not a gas phase
    if isnan(v1+v2)
        return v1 + v2
    end
    #the test was done on equal volume chambers (300 cc), but mathematically it doesn't matter
    v_test = 1.0
    z1 = v_test/v1
    z2 = v_test/v2
    x = [z1,z2]
    x ./= sum(x)
    B̄ = B(model,T,x)
    
    B12 = (B̄ - x[1]*x[1]*B11 - x[2]*x[2]*B22)/(2*x[1]*x[2])    
    return B12
end


function VT_compressibility_factor(model::EoSModel, V, T, z=SA[1.])
    p = pressure(model,V,T,z)
    return p*V/(sum(z)*R̄*T)
end

"""
    pip(model::EoSModel,V,T,z=[1.0])

Phase identification parameter `Π`. as described in _1_. If `Π > 1`, then the phase is clasified as a liquid or a liquid-like vapor, being a vapor or vapor-like liquid otherwise.

This identification parameter fails at temperatures and pressures well aboVe the critical point.

Calculated as:
```
Π = V*((∂²p/∂V∂T)/(∂p/∂T) - (∂²p/∂V²)/(∂p/∂V))
```


1.  G. Venkatarathnama, L.R. Oellrich, Identification of the phase of a fluid using partial derivatives of pressure, volume,and temperature without reference to saturation properties: Applications in phase equilibria calculations, Fluid Phase Equilibria 301 (2011) 225–233


"""
function pip(model::EoSModel, V, T, z=SA[1.0])
    _∂2p = ∂2p(model,V,T,z)
    hess_p, grad_p, _ = _∂2p
    Π = V*(hess_p[1,2]/grad_p[2]  - hess_p[1,1]/grad_p[1])
end

function VT_mass_density(model::EoSModel,V,T,z=SA[1.0])
    molar_weight = molecular_weight(model,z)
    return molar_weight/V
end

function VT_molar_density(model::EoSModel,V,T,z=SA[1.0])
    return sum(z)/V
end

#Vector Properties

function VT_partial_property(model::EoSModel,V,T,z,property::ℜ) where {ℜ}
    fun(x) = property(model,V,T,x)
    TT = gradient_type(V,T,z)
    return ForwardDiff.gradient(fun,z)::TT
end

VT_chemical_potential(model::EoSModel, V, T, z=SA[1.]) = VT_partial_property(model,V,T,z,eos)
VT_chemical_potential_res(model::EoSModel, V, T, z=SA[1.]) = VT_partial_property(model,V,T,z,eos_res)

export second_virial_coefficient,pressure,cross_second_virial,equivol_cross_second_virial

