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

function VT_chemical_potential(model::EoSModel, V, T, z=SA[1.])
    fun(x) = eos(model,V,T,x)
    return ForwardDiff.gradient(fun,z)
end

function VT_chemical_potential_res(model::EoSModel, V, T, z=SA[1.])
    fun(x) = eos_res(model,V,T,x)
    return ForwardDiff.gradient(fun,z)
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

function VT_helmholtz_free_energy(model::EoSModel, V, T, z=SA[1.])
    return eos(model,V,T,z)
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
    p0,dpdV = p∂p∂V(model,V,T,z)
    return -1/V/dpdV
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
    return ∂²A∂V∂T/(V*∂²A∂V²)
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

Calculates the second virial coefficient `B`, defined as:

```julia
B = lim(V->∞)[ V^2/RT *  (∂Aᵣ∂V + V*∂²Aᵣ∂V²) ]
```
where `Aᵣ` is the residual helmholtz energy.
"""
function second_virial_coefficient(model::EoSModel, T, z=SA[1.])
    TT = promote_type(eltype(z),typeof(T))
    V = 1/sqrt(eps(TT))
    fAᵣ(x) = eos_res(model,x,T,z)
    Aᵣ,∂Aᵣ∂V,∂²Aᵣ∂V² = Solvers.f∂f∂2f(fAᵣ,V)
    return V^2/(R̄*T)*(∂Aᵣ∂V+V*∂²Aᵣ∂V²)
end

function VT_compressibility_factor(model::EoSModel, V, T, z=SA[1.])
    p = pressure(model,V,T,z)
    return p*V/(R̄*T)
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

export second_virial_coefficient,pressure
