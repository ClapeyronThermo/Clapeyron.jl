function pressure(model::EoSModel, V, T, z=SA[1.])
    return -∂f∂V(model,V,T,z)
end

function VT_entropy(model::EoSModel, V, T, z=SA[1.])
    return -∂f∂T(model,V,T,z)
end

function VT_chemical_potential(model::EoSModel, V, T, z=SA[1.])
    fun(x) = eos(model,V,T,x)
    return ForwardDiff.gradient(fun,z)
end

function VT_internal_energy(model::EoSModel, V, T, z=SA[1.])
    _df,f = ∂f(model,V,T,z)
    dV,dT = _df
    return f - dT*T
end

function VT_enthalpy(model::EoSModel, V, T, z=SA[1.])
    _df,f = ∂f(model,V,T,z)
    dV,dT = _df
    return f - dV*V - dT*T
end

function VT_gibbs_free_energy(model::EoSModel, V, T, z=SA[1.])
    _df,f = ∂f(model,V,T,z)
    dV,dT = _df
    return f - dV*V
end

function VT_helmholtz_free_energy(model::EoSModel, V, T, z=SA[1.])
    return eos(model,V,T,z)
end

function VT_isochoric_heat_capacity(model::EoSModel, V, T, z=SA[1.])
    fun(x) = eos(model,V,x,z)
    df(x) = ForwardDiff.derivative(fun,x)
    d2f(x) = ForwardDiff.derivative(df,x)
    return -T*d2f(T)
end

function VT_isobaric_heat_capacity(model::EoSModel, V, T, z=SA[1.])
    d2f = f_hess(model,V,T,z)
    return T*(d2f[1,2]^2/d2f[1]-d2f[2,2])
end

function VT_isothermal_compressibility(model::EoSModel, V, T, z=SA[1.])
    p0,dpdV = p∂p∂V(model,V,T,z)
    return -1/V*dpdV^-1
end

function VT_isentropic_compressibility(model::EoSModel, V, T, z=SA[1.])
    d2f = f_hess(model,V,T,z)
    return 1/V*(d2f[1]-d2f[1,2]^2/d2f[2,2])^-1
end

function VT_speed_of_sound(model::EoSModel, V, T, z=SA[1.])
    Mr = molecular_weight(model,z)
    d2f = f_hess(model,V,T,z)
    return V*sqrt((d2f[1]-d2f[1,2]^2/d2f[2,2])/Mr)
end

function VT_isobaric_expansivity(model::EoSModel, V, T, z=SA[1.])
    d2f = f_hess(model,V,T,z)
    return d2f[1,2]/(V*d2f[1])
end

function VT_joule_thomson_coefficient(model::EoSModel, V, T, z=SA[1.])
    d2f = f_hess(model,V,T,z)
    return -(d2f[1,2]-d2f[1]*((T*d2f[2,2]+V*d2f[1,2])/(T*d2f[1,2]+V*d2f[1])))^-1
end

function second_virial_coefficient(model::EoSModel, T, z=SA[1.])
    TT = promote_type(eltype(z),typeof(T))
    V = 1/sqrt(eps(TT))
    fun(x) = eos_res(model,x[1],T,z)
    df(x) = ForwardDiff.derivative(fun,x[1])
    d2f(x) = ForwardDiff.derivative(df,x[1])
    return V^2/(R̄*T)*(df(V)+V*d2f(V))
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

export second_virial_coefficient
