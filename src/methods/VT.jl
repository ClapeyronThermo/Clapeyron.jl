"""
    pressure(model::EoSModel, V, T, z=SA[1.])

Default units: `[Pa]`

Returns the pressure of the model at a given volume `V`, temperature `T` and composition `z`, defined as:

```julia
p = -∂A/∂V
```
where A is the Helmholtz energy `[J]`,
V is the volume `[m³]`
"""
function pressure(model::EoSModel, V, T, z=SA[1.])
    return VT_pressure(model, V, T, z)
end

function temperature end

VT_pressure(model, V, T) = VT_pressure(model, V, T, SA[1.0])
VT_pressure(model, V, T, z) = -∂f∂V(model,V,T,z)
VT_temperature(model, V, T, z=SA[1.]) = T
VT_volume(model, V, T, z=SA[1.]) = V

function pressure_res(model::EoSModel, V, T, z=SA[1.])
    fun(x) = eos_res(model,x,T,z)
    return -Solvers.derivative(fun,V)
end
VT_pressure_res(model, V, T) = VT_pressure_res(model, V, T, SA[1.])
VT_pressure_res(model, V, T, z) = pressure_res(model,V,T,z)

function VT_entropy(model::EoSModel, V, T, z::AbstractVector=SA[1.])
    return -∂f∂T(model,V,T,z)
end

VT_mass_entropy(model::EoSModel,V, T, z::AbstractVector = SA[1.0]) = VT_entropy(model,V,T,z)/molecular_weight(model,z)


function VT_entropy_res(model::EoSModel, V, T, z=SA[1.])
    fun(x) = eos_res(model,V,x,z)
    return -Solvers.derivative(fun,T)
end

function VT_internal_energy(model::EoSModel, V, T, z::AbstractVector=SA[1.])
    ideal = model isa IdealModel
    if iszero(1/V) && !ideal
        return VT_internal_energy(idealmodel(model),V,T,z)
    end

    if ideal
        V₀ = oneunit(V) #the volume term gets cancelled out
    else
        V₀ = V
    end

    A, ∂A∂T = f∂fdT(model,V₀,T,z)
    return A - T*∂A∂T
end

VT_mass_internal_energy(model::EoSModel,V, T, z::AbstractVector = SA[1.0]) = VT_internal_energy(model,V,T,z)/molecular_weight(model,z)


function VT_internal_energy_res(model::EoSModel, V, T, z=SA[1.])
    A, ∂A∂V, ∂A∂T = ∂f_res_vec(model,V,T,z)
    return A - T*∂A∂T
end

function VT_enthalpy(model::EoSModel, V, T, z::AbstractVector=SA[1.])
    ideal = model isa IdealModel
    if iszero(1/V) && !ideal
        return VT_internal_energy(idealmodel(model),V,T,z)
    end

    if model isa IdealModel
        V₀ = oneunit(V) #the volume term gets cancelled out
    else
        V₀ = V
    end

    A, ∂A∂V, ∂A∂T = ∂f_vec(model,V₀,T,z)

    if ideal
        return A + sum(z)*Rgas(model)*T - T*∂A∂T
    else
        return A - V*∂A∂V - T*∂A∂T
    end
end

VT_mass_enthalpy(model::EoSModel,V, T, z::AbstractVector = SA[1.0]) = VT_enthalpy(model,V,T,z)/molecular_weight(model,z)

function VT_enthalpy_res(model::EoSModel, V, T, z=SA[1.])
    A, ∂A∂V, ∂A∂T = ∂f_res_vec(model,V,T,z)
    PrV = ifelse(isinf(primalval(V)),zero(∂A∂V),- V*∂A∂V)
    return A + PrV - T*∂A∂T
end

function VT_gibbs_free_energy(model::EoSModel, V, T, z::AbstractVector=SA[1.], p = nothing)
    ideal = model isa IdealModel
    if iszero(1/V) && !ideal
        return VT_gibbs_free_energy(idealmodel(model), V, T, z)
    end

    if p == nothing
        A,∂A∂V = f∂fdV(model,V,T,z)
    else
        A = eos(model,V,T,z)
        ∂A∂V = -p
    end

    if ideal
        return A + sum(z)*Rgas(model)*T
    else
        return A - V*∂A∂V
    end
end

VT_mass_gibbs_free_energy(model::EoSModel,V, T, z::AbstractVector = SA[1.0],p = nothing) = VT_gibbs_free_energy(model,V,T,z,p)/molecular_weight(model,z)

function VT_gibbs_free_energy_res(model::EoSModel, V, T, z=SA[1.])
    fun(x) = eos_res(model,x,T,z)
    Ar,∂A∂Vr = Solvers.f∂f(fun,V)
    PrV = ifelse(iszero(1/V),zero(∂A∂Vr),- V*∂A∂Vr)
    return Ar + PrV
end

function VT_helmholtz_free_energy(model::EoSModel, V, T, z::AbstractVector=SA[1.])
    return eos(model,V,T,z)
end

VT_mass_helmholtz_free_energy(model::EoSModel,V, T, z::AbstractVector = SA[1.0]) = VT_helmholtz_free_energy(model,V,T,z)/molecular_weight(model,z)


function VT_helmholtz_free_energy_res(model::EoSModel, V, T, z=SA[1.])
    return eos_res(model,V,T,z)
end

const VT_helmholtz_energy = VT_helmholtz_free_energy
const VT_helmholtz_energy_res = VT_helmholtz_free_energy_res
const VT_gibbs_energy = VT_gibbs_free_energy
const VT_gibbs_energy_res = VT_gibbs_free_energy_res

function VT_isochoric_heat_capacity(model::EoSModel, V, T, z=SA[1.])
    ∂²A∂T² = ∂²f∂T²(model,V,T,z)
    return -T*∂²A∂T²
end

VT_mass_isochoric_heat_capacity(model::EoSModel,V, T, z::AbstractVector = SA[1.0]) = VT_isochoric_heat_capacity(model,V,T,z)/molecular_weight(model,z)


function VT_isobaric_heat_capacity(model::EoSModel, V, T, z=SA[1.])
    if iszero(1/V) || model isa IdealModel
        ∂²A∂T² = ∂²f∂T²(model,V,T,z)
        return -T*∂²A∂T² + Rgas(model)*sum(z)
    else
        d²A = f_hess(model,V,T,z)
        ∂²A∂V∂T = d²A[1,2]
        ∂²A∂V² = d²A[1,1]
        ∂²A∂T² = d²A[2,2]
        return -T*(∂²A∂T² - ∂²A∂V∂T^2/∂²A∂V²)
    end
end

VT_mass_isobaric_heat_capacity(model::EoSModel,V, T, z::AbstractVector = SA[1.0]) = VT_isobaric_heat_capacity(model,V,T,z)/molecular_weight(model,z)

function VT_adiabatic_index(model::EoSModel, V, T, z=SA[1.])
    if iszero(1/V) || model isa IdealModel
        ∂²A∂T² = ∂²f∂T²(model,V,T,z)
        1 - Rgas(model)*sum(z)/(∂²A∂T²*T)
    else
        d²A = f_hess(model,V,T,z)
        ∂²A∂V∂T = d²A[1,2]
        ∂²A∂V² = d²A[1,1]
        ∂²A∂T² = d²A[2,2]
        return 1 - ∂²A∂V∂T*∂²A∂V∂T/(∂²A∂V²*∂²A∂T²)
    end
end

function VT_isothermal_compressibility(model::EoSModel, V, T, z=SA[1.])
    if iszero(1/V) || model isa IdealModel
        return V/(sum(z)*Rgas(model)*T)
    else
        _,∂p∂V = p∂p∂V(model,V,T,z)
        return -1/V/∂p∂V
    end
end

function VT_isentropic_compressibility(model::EoSModel, V, T, z=SA[1.])
    if iszero(1/V) || model isa IdealModel
        ∂²A∂T² = ∂²f∂T²(model,V,T,z)
        R = Rgas(model)
        V_∂²A∂V∂T_2 = R*R/V
        V_∂²A∂V² = R*T*sum(z)/V
        return 1/(V_∂²A∂V² - V_∂²A∂V∂T_2/∂²A∂T²)
    else
        d²A = f_hess(model,V,T,z)
        ∂²A∂V∂T = d²A[1,2]
        ∂²A∂V² = d²A[1,1]
        ∂²A∂T² = d²A[2,2]
        return 1/V/(∂²A∂V²-∂²A∂V∂T^2/∂²A∂T²)
    end
end

function VT_speed_of_sound(model::EoSModel, V, T, z=SA[1.])
    Mr = molecular_weight(model,z)
    if iszero(1/V) || model isa IdealModel
        γ = VT_adiabatic_index(model,V,T,z)
        return sqrt(γ*Rgas(model)*T*sum(z)/Mr)
    else
        d²A = f_hess(model,V,T,z)
        ∂²A∂V∂T = d²A[1,2]
        ∂²A∂V² = d²A[1,1]
        ∂²A∂T² = d²A[2,2]
        return V*sqrt((∂²A∂V²-∂²A∂V∂T^2/∂²A∂T²)/Mr)
    end
end

function VT_isobaric_expansivity(model::EoSModel, V, T, z=SA[1.])
    if iszero(1/V)  || model isa IdealModel
        return one(Base.promote_eltype(model,V,T,z))/T
    else
        d²A = f_hess(model,V,T,z)
        ∂²A∂V∂T = d²A[1,2]
        ∂²A∂V² = d²A[1,1]
        return -∂²A∂V∂T/(V*∂²A∂V²)
    end
end

function VT_joule_thomson_coefficient(model::EoSModel, V, T, z=SA[1.])
    if iszero(1/V) || model isa IdealModel
        B,∂B∂T = B∂B∂T(model,T,z)
        Cp = VT_isobaric_heat_capacity(model,V,T,z)
        return (T*∂B∂T - B)/Cp
    else
        d²A = f_hess(model,V,T,z)
        ∂²A∂V∂T = d²A[1,2]
        ∂²A∂V² = d²A[1,1]
        ∂²A∂T² = d²A[2,2]
        return -(∂²A∂V∂T - ∂²A∂V²*((T*∂²A∂T² + V*∂²A∂V∂T) / (T*∂²A∂V∂T + V*∂²A∂V²)))^-1
    end
end

function VT_compressibility_factor(model::EoSModel, V, T, z=SA[1.],p = nothing)
    if p === nothing
        return pressure(model,V,T,z)*V/(sum(z)*R̄*T)
    else
        return p*V/(sum(z)*R̄*T)
    end
end

"""
    second_virial_coefficient(model::EoSModel, T, z=SA[1.])

Default units: `[m³]`

Calculates the second virial coefficient `B`, defined as:

```julia
B = lim(ρ->0)[∂Aᵣ/∂ρ]
```
where `Aᵣ` is the residual Helmholtz energy.
"""
function second_virial_coefficient(model::EoSModel, T, z=SA[1.])
   return second_virial_coefficient_impl(model,T,z)
end

function second_virial_coefficient_impl(model::EoSModel, T, z = SA[1.0])
    TT = one(Base.promote_eltype(model,T,z))
    ϵ = 1/eps(TT)
    V = sqrt(sum(z)*ϵ)
    return pressure_res(model,V,T,z)*ϵ/(Rgas(model)*T)
end

function B∂B∂T(model,T,z = SA[1.0])
    b(T) = second_virial_coefficient(model,T,z)
    return Solvers.f∂f(b,T)
end
"""
    cross_second_virial(model,T,z)

Default units: `[m³]`

Calculates the second cross virial coefficient (B₁₂) of a binary mixture, using the definition:

```julia
B̄ = x₁^2*B₁₁ + 2x₁x₂B₁₂ + x₂^2*B₂₂
B₁₂ = (B̄ - x₁^2*B₁₁ - x₂^2*B₂₂)/2x₁x₂
```


!!! info "Composition-dependent property"
    The second cross virial coefficient calculated from an equation of state can present a dependency on composition [1], but normally, experiments for obtaining the second virial coefficient are made by mixing the same volume of two gases. You can calculate B₁₂ in this way by using (Clapeyron.equivol_cross_second_virial)[@ref]

## References
1. Jäger, A., Breitkopf, C., & Richter, M. (2021). The representation of cross second virial coefficients by multifluid mixture models and other equations of state. Industrial & Engineering Chemistry Research, 60(25), 9286–9295. [doi:10.1021/acs.iecr.1c01186](https://doi.org/10.1021/acs.iecr.1c01186)
"""
function cross_second_virial(model,T,z)
    B = second_virial_coefficient
    n = length(model)
    ∑z = sum(z)
    if n == 1
        return zero(T + first(z))
    else
        binary_component_check(cross_second_virial,model)
        model1,model2 = split_pure_model(model)
        B̄ = B(model,T,z)/∑z #1 mol
        B1,B2 = B(model1,T),B(model2,T) #1 mol by default
        #@show B1,B2,B̄
        x = z/∑z
        #B̄ = (B1*x1^2 + B2*x2^2 + B12*x1*x2)
        B12 = (B̄ - x[1]*x[1]*B1 - x[2]*x[2]*B2)/(2*x[1]*x[2])
        return B12*∑z
    end
end

"""
    equivol_cross_second_virial(model::EoSModel,T,p_exp = 200000.0)

Calculates the second cross virial coefficient, by simulating the mixing of equal volumes of pure gas, at T,P conditions.
The equal volume of each pure gas sets an specific molar amount for each component. Details of the experiment can be found at [1].

## Example
```
model = SAFTVRQMie(["helium","neon"])
B12 = equivol_cross_second_virial(model,model,T,p_exp = 200000.0)

```
## References
1. Brewer, J., & Vaughn, G. W. (1969). Measurement and correlation of some interaction second virial coefficients from − 125° to 50°C. I. The Journal of Chemical Physics, 50(7), 2960–2968. [doi:10.1063/1.1671491](https://doi.org/10.1063/1.1671491)
"""
function equivol_cross_second_virial(model,T,p_exp = 200000.0)
    binary_component_check(cross_second_virial,model)
    #they do experiments at constant volume and temperature, so we are gonna need to calculate the mole fractions for that
    m1,m2 = split_pure_model(model)
    B = second_virial_coefficient
    B11,B22 = B(m1,T),B(m2,T)
    v01,v02 = volume_virial(B11,p_exp,T),volume_virial(B22,p_exp,T)
    v1,v2 = volume(m1,p_exp,T,vol0 = v01),volume(m1,p_exp,T,vol0 = v02) #we do this in case there is not a gas phase
    if isnan(v1+v2)
        return v1 + v2
    end
    #the test was done on equal volume chambers (300 cm³), but mathematically it doesn't matter
    v_test = 1.0
    z1 = v_test/v1
    z2 = v_test/v2
    x = [z1,z2]
    x ./= sum(x)
    B̄ = B(model,T,x)

    B12 = (B̄ - x[1]*x[1]*B11 - x[2]*x[2]*B22)/(2*x[1]*x[2])
    return B12
end

"""
    pip(model::EoSModel,V,T,z=[1.0])

Phase identification parameter `Π`, as described in _1_. If `Π > 1`, then the phase is clasified as a liquid or a liquid-like vapor, being a vapor or vapor-like liquid otherwise.

This identification parameter fails at temperatures and pressures well above the critical point.

Calculated as:
```
Π = V*((∂²p/∂V∂T)/(∂p/∂T) - (∂²p/∂V²)/(∂p/∂V))
```
## References
1.  G. Venkatarathnama, L.R. Oellrich, Identification of the phase of a fluid using partial derivatives of pressure, volume,and temperature without reference to saturation properties: Applications in phase equilibria calculations, Fluid Phase Equilibria 301 (2011) 225–233
"""
function pip(model::EoSModel, V, T, z=SA[1.0])
    Π,∂p∂V = _pip(model,V,T,z)
    return Π
end

VT_pip(model::EoSModel, V, T, z=SA[1.0]) = pip(model,V,T,z)

function _pip(model::EoSModel, V, T, z=SA[1.0])
    _∂2p = ∂2p(model,V,T,z)
    hess_p, grad_p, _ = _∂2p
    ∂p∂V = grad_p[1]
    Π = V*(hess_p[1,2]/grad_p[2]  - hess_p[1,1]/grad_p[1])
    return Π,∂p∂V
end

function VT_fundamental_derivative_of_gas_dynamics(model::EoSModel, V, T, z=SA[1.0])
    Mr = molecular_weight(model,z)
    d²A = f_hess(model,V,T,z)
    ∂²A∂V∂T = d²A[1,2]
    ∂²A∂V² = d²A[1,1]
    ∂²A∂T² = d²A[2,2]
    c =  V*sqrt((∂²A∂V²-∂²A∂V∂T^2/∂²A∂T²)/Mr)
    A(x) = eos(model,V,x,z)
    ∂A∂T(x) = Solvers.derivative(A,x)
    ∂²A∂T²(x) = -T*Solvers.derivative(∂A∂T,x)
    Cᵥ,∂Cᵥ∂T = Solvers.f∂f(∂²A∂T²,T)
    _∂2p = ∂2p(model,V,T,z)
    hess_p, grad_p, _ = _∂2p
    ∂²p∂T²,∂²p∂V²,∂²p∂V∂T = hess_p[2,2],hess_p[1,1],hess_p[1,2]
    ∂p∂T,∂p∂V = grad_p[2],grad_p[1]
    Γ₁ = ∂²p∂V²
    Γ₂ = (-3*T/Cᵥ)*∂p∂T*∂²p∂V∂T
    Γ₃ = ((T/Cᵥ)*∂p∂T)^2 * (3*∂²p∂T² + (∂p∂T/T)*(1 - (T/Cᵥ)*∂Cᵥ∂T))
    return (V*V*V/(2*c*c*Mr))*(Γ₁ + Γ₂ + Γ₃)
end

"""
    VT_identify_phase(model::EoSModel, V, T, z=SA[1.0])::Symbol

Returns the phase of a fluid at the conditions specified by volume `V`, temperature`T` and composition `z`.
Uses the phase identification parameter criteria from `Clapeyron.pip`

Returns `liquid` if the phase is liquid (or liquid-like), `vapour` if the phase is vapour (or vapour-like), and `:unknown` if the calculation of the phase identification parameter failed (the V-T-z point was mechanically unstable).
"""
function VT_identify_phase(model::EoSModel, V, T, z=SA[1.0])
    Π,∂p∂V = _pip(model, V, T, z)
    βT = -1/V/∂p∂V
    if Π > 1 && βT >= 0
        return :liquid
    elseif Π <= 1 && βT >= 0
        return :vapour
    else #the calculation failed
        return :unknown
    end
end

function VT_mass_density(model::EoSModel,V,T,z=SA[1.0])
    molar_weight = molecular_weight(model,z)
    return molar_weight/V
end

function VT_molar_density(model::EoSModel,V,T,z=SA[1.0])
    return sum(z)/V
end

#Vector Properties

struct ZVar{P,M,V,T}
    property::P
    model::M
    vol::V
    temp::T
end

(fixed::ZVar{P,M,V,T})(z::Z) where {P,M,V,T,Z} = fixed.property(fixed.model,fixed.vol,fixed.temp,z)

function VT_molar_gradient(model::EoSModel,V,T,z::AbstractVector,property::ℜ) where {ℜ}
    fun = ZVar(property,model,V,T)
    TT = gradient_type(model,T+V,z)
    return Solvers.gradient(fun,z)::TT
end

function VT_molar_gradient!(fx::F,model::EoSModel,V,T,z,property::ℜ) where {F<:AbstractVector,ℜ}
    if isnan(V) || isnan(T) || any(isnan,z)
        fx .= NaN
        return fx
    end
    fun = ZVar(property,model,V,T)
    return Solvers.gradient!(fx,fun,z)::F
end

function VT_molar_gradient!(cache::F,model::EoSModel,V,T,z,property::ℜ) where {F<:Tuple,ℜ}

    result,aux,∇f,A1,x1,x2,x3,hconfig = cache
    if isnan(V) || isnan(T) || any(isnan,z)
        ∇f .= NaN
        return ∇f
    end
    nc = length(z)
    if nc == 1
        f1(_z) = property(model,V,T,SVector(_z))
        ∇f[1] = ForwardDiff.derivative(f1,z[1])
    else
        aux .= 0
        aux[1:nc] = z
        gconfig = Solvers._GradientConfig(hconfig)
        fun(aux) = property(model, V, T, @view(aux[1:nc]))
        _result = ForwardDiff.gradient!(result, fun, aux, gconfig, Val{false}())
        dresult = DiffResults.gradient(_result)
        fx_temp = @view dresult[1:nc]
        ∇f .= fx_temp
    end
    return ∇f
end

VT_chemical_potential(model::EoSModel, V, T, z=SA[1.]) = VT_molar_gradient(model,V,T,z,eos)
VT_chemical_potential_res(model::EoSModel, V, T, z=SA[1.]) = VT_molar_gradient(model,V,T,z,eos_res)
VT_chemical_potential_res!(r, model::EoSModel, V, T, z=SA[1.]) = VT_molar_gradient!(r,model,V,T,z,eos_res)
VT_chemical_potential!(result,model,V,T,z) = VT_molar_gradient!(result,model,V,T,z,eos)

function VT_fugacity_coefficient(model::EoSModel,V,T,z=SA[1.])
    return _VT_fugacity_coefficient(model,V,T,z)
end

#i saw some code that calls only(fugacity_coefficient(model,p,T)).
#this is specialization for the case of a single component.
function _VT_fugacity_coefficient(model::EoSModel,V,T,z)
    p = pressure(model,V,T,z)
    μ_res = VT_chemical_potential_res(model,V,T,z)
    R̄ = Rgas(model)
    Z = p*V/R̄/T/sum(z)
    return exp.(μ_res ./ R̄ ./ T) ./ Z
end

function _VT_fugacity_coefficient(model::EoSModel,V,T,z::SingleComp)
    f(_V) = eos_res(model, _V, T,z)
    A,dAdV = Solvers.f∂f(f,V)
    R̄ = Rgas(model)
    ∑z= sum(z)
    p_ideal = ∑z*R̄*T/V
    p = -dAdV + p_ideal
    μ_res = muladd(-V,dAdV,A)
    Z = p*V/R̄/T/sum(z)
    ϕ = exp(μ_res/R̄/T)/Z
    return SVector(ϕ)
end

function VT_fugacity_coefficient!(φ,model::EoSModel,V,T,z=SA[1.],p = pressure(model,V,T,z))
    if isnan(V)
        φ .= NaN
        return φ
    end
    φ = VT_chemical_potential_res!(φ,model,V,T,z)
    R̄ = Rgas(model)
    Z = p*V/R̄/T/sum(z)
    φ ./= (R̄*T)
    φ .= exp.(φ)
    φ ./= Z
    return φ
end

function VT_thermodynamic_factor(model::EoSModel, V, T, z)
    N = length(model)
    ∑z = sum(z)
    v = V / ∑z
    x = z ./ ∑z
    xN1 = @view x[1:N-1]
    RT = Rgas(model) * T

    fun_A(_xv) = begin
        _x = @view _xv[1:N]
        _v = last(_xv)
        return eos(model, _v, T, _x)
    end
    H_A = Solvers.hessian(fun_A, vcat(x, v))
    H_A_nn = @view H_A[1:N,1:N]
    H_A_nv = @view H_A[N+1,1:N] 
    H_A_vv⁻¹ = inv(H_A[N+1,N+1])
    H_G_nn = H_A_nn .- (H_A_nv * H_A_nv') .* H_A_vv⁻¹

    H_G_nN = @view H_G_nn[1:N-1,N]
    H_G_Nn = @view H_G_nn[N,1:N-1]
    H_G_xx = H_G_nn[1:N-1,1:N-1] .- H_G_nN .- H_G_Nn .+ H_G_nn[N,N]
    F = [i == j ? 1-x[i] : -x[i] for i in 1:N-1, j in 1:N-1]
    ∂μᵢ∂xⱼ = H_G_xx * F

    Γ = xN1 ./ RT .* ∂μᵢ∂xⱼ
    return Γ
end

const VT_helmholtz_energy = VT_helmholtz_free_energy
const VT_gibbs_energy = VT_gibbs_free_energy
const VT_mass_helmholtz_energy = VT_mass_helmholtz_free_energy
const VT_mass_gibbs_energy = VT_mass_gibbs_free_energy


export pressure
export second_virial_coefficient,cross_second_virial,equivol_cross_second_virial

const CLAPEYRON_PROPS = [:temperature,:volume, :pressure, :entropy, :internal_energy, :enthalpy, :gibbs_free_energy, :helmholtz_free_energy,
                    :entropy_res, :internal_energy_res, :enthalpy_res, :gibbs_free_energy_res, :helmholtz_free_energy_res,
                    :helmholtz_energy,:gibbs_energy,
                    #mass properties, first order
                    :mass_entropy,:mass_enthalpy,:mass_internal_energy,:mass_gibbs_free_energy,:mass_helmholtz_free_energy,
                    :mass_helmholtz_energy,:mass_gibbs_energy,
                    #second derivative order properties
                    :isochoric_heat_capacity, :isobaric_heat_capacity, :adiabatic_index,
                    :isothermal_compressibility, :isentropic_compressibility, :speed_of_sound,
                    :isobaric_expansivity, :joule_thomson_coefficient, :inversion_temperature,
                    #second derivative order, mass properties
                    :mass_isobaric_heat_capacity,:mass_isochoric_heat_capacity,
                    #higher derivative order properties
                    :fundamental_derivative_of_gas_dynamics,
                    #volume properties
                    :mass_density, :molar_density, :compressibility_factor,
                    #other properties
                    :identify_phase, :pip,
]

const CLAPEYRON_PROP_ALIASES = [:mass_gibbs_energy,:gibbs_energy,:helmholtz_energy,:mass_helmholtz_energy]

function VT_symbol(x::Symbol)
    return Symbol(:VT_,x)
end


#module used to translate between the normal symbol and the VT_symbol.
module VT0
    using Clapeyron: Clapeyron, CLAPEYRON_PROPS
    for prop in Clapeyron.CLAPEYRON_PROPS
        VT_prop = Clapeyron.VT_symbol(prop)
        @eval begin
            function $prop(model,V,T,z=Clapeyron.SA[1.])
                return Clapeyron.$VT_prop(model,V,T,z)
            end
        end
    end
    chemical_potential(model,V,T,z=Clapeyron.SA[1.0]) = Clapeyron.VT_chemical_potential(model,V,T,z)
    chemical_potential_res(model,V,T,z = Clapeyron.SA[1.0]) = Clapeyron.VT_chemical_potential_res(model,V,T,z)
end
