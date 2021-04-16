function pressure(model::EoSModel, v, T,  z=SA[1.])
    return -∂f∂v(model,v,T,z)
end

function vt_entropy(model::EoSModel, v, T,   z=SA[1.])
    return -∂f∂t(model,v,T,z)
end

function vt_chemical_potential(model::EoSModel, v, T, z= SA[1.])
    fun(x) = eos(model, v, T,z)
    return ForwardDiff.gradient(fun,z)
end

function vt_internal_energy(model::EoSModel, v, T,  z=SA[1.])
    _df,f =  ∂f(model,v,T,z)
    dv,dt = _df
    return f  - dt*T
end

function vt_enthalpy(model::EoSModel, v, T,  z=SA[1.])
    _df,f =  ∂f(model,v,T,z)
    dv,dt = _df
    return f  - dv*v - dt*T
end

function vt_gibbs_free_energy(model::EoSModel, v, T,  z=SA[1.])
    _df,f =  ∂f(model,v,T,z)
    dv,dt = _df
    return f  - dv*v
end

function vt_helmholtz_free_energy(model::EoSModel, v, T,  z=SA[1.])
    return eos(model, v, T,z)
end

function vt_isochoric_heat_capacity(model::EoSModel, v, T,  z=SA[1.])
    fun(x)  = eos(model, v, x,z)
    df(x)   = ForwardDiff.derivative(fun,x)
    d2f(x)  = ForwardDiff.derivative(df,x)
    return -T*d2f(T)
end

function vt_isobaric_heat_capacity(model::EoSModel, v, T,  z=SA[1.])
    d2f = f_hess(model,v,T,z)
    return T*(d2f[1,2]^2/d2f[1]-d2f[2,2])
end

function vt_isothermal_compressibility(model::EoSModel, v, T,  z=SA[1.])
    p0,dpdv = p∂p∂v(model,v,T,z)
    return -1/v*dpdv^-1
end

function vt_isentropic_compressibility(model::EoSModel, v, T,  z=SA[1.])
    d2f = f_hess(model,v,T,z)
    return 1/v*(d2f[1]-d2f[1,2]^2/d2f[2,2])^-1
end

function vt_speed_of_sound(model::EoSModel, v, T,  z=SA[1.])
    Mr      = molecular_weight(model,z)
    d2f = f_hess(model,v,T,z)
    return v*sqrt((d2f[1]-d2f[1,2]^2/d2f[2,2])/Mr)
end

function vt_isobaric_expansivity(model::EoSModel, v, T,  z=SA[1.])
    d2f = f_hess(model,v,T,z)
    return d2f[1,2]/(v*d2f[1])
end

function vt_joule_thomson_coefficient(model::EoSModel, v, T,  z=SA[1.])

    d2f = f_hess(model,v,T,z)
    return -(d2f[1,2]-d2f[1]*((T*d2f[2,2]+v*d2f[1,2])/(T*d2f[1,2]+v*d2f[1])))^-1
end


function second_virial_coefficient(model::EoSModel, T,  z=SA[1.])
    TT = promote_type(eltype(z),typeof(T))  
    V = 1/sqrt(eps(TT))
    _∂2f = ∂2f(model,V,T,z)
    hessf,gradf,f = _∂2f
    _p,dpdv = p∂p∂v(model,V,T,z)
    return -V^2/(R̄*T)*(_p+V*dpdv)
end

function vt_compressibility_factor(model::EoSModel, v, T,  z=SA[1.])
    p = pressure(model,v,T,z)
    return p*v/(R̄*T)
end

"""
    pip(model::EoSModel,v,T,z=[1.0])

Phase identification parameter `Π`. as described in _1_. If `Π > 1`, then the phase is clasified as a liquid or a liquid-like vapor, being a vapor or vapor-like liquid otherwise. 

This identification parameter fails at temperatures and pressures well above the critical point.

Calculated as:
```
Π = v*((∂²p/∂v∂T)/(∂p/∂T) - (∂²p/∂v²)/(∂p/∂v)) 
```


1.  G. Venkatarathnama, L.R. Oellrich, Identification of the phase of a fluid using partial derivatives of pressure, volume,and temperature without reference to saturation properties: Applications in phase equilibria calculations, Fluid Phase Equilibria 301 (2011) 225–233


"""
function pip(model::EoSModel,v,T,z=SA[1.0])
    _∂2p = ∂2p(model,v,T,z)
    hess_p, grad_p, _ = _∂2p
    Π = v*(hess_p[1,2]/grad_p[2]  - hess_p[1,1]/grad_p[1]) 
end
