struct Ice2009 <: GibbsBasedModel
    components::Vector{String}
    s0::Float64
    Rgas::Float64
    references::Vector{String}

end

"""
    Ice2009 <: EmpiricHelmholtzModel
    Ice2009()

## Input parameters

None

## Description

IAPWS (International Association for the Properties of Water and Steam) Pure water Model for ice Ih, 2009 update.

## References

1. Feistel, R., & Wagner, W. (2006). A new equation of state for H2O ice Ih. Journal of Physical and Chemical Reference Data, 35(2), 1021–1047. [doi:10.1063/1.2183324](https://doi.org/10.1063/1.2183324)
2. IAPWS R10-06 (2009). Revised Release on the Equation of State 2006 for H2O Ice Ih
"""
Ice2009

function Ice2009(components = ["water"];s0 = -0.332733756492168e4,Rgas = 8.314472, idealmodel = nothing,verbose = false,reference_state = nothing)
    _components = format_components(components)
    references = default_references(Ice2009)
    return Ice2009(_components,s0,Rgas,references)
end

Rgas(model::Ice2009) = model.Rgas
default_references(::Type{Ice2009}) = ["IAPWS R10-06(2009)","10.1063/1.2183324"]
is_splittable(model::Ice2009) = false

function eos_g(model::Ice2009,p,T,z)
    #Δπ = (p - 101325.0)/611.654771007894
    Δπ = (p - 101325.0)/611.657
    Tt = 273.16
    τ = T/Tt
    s0 = model.s0

    r1 = 0.447050716285388e2 + 0.656876847463481e2im
    r2 = ice_r2(model,Δπ)

    t1 = 0.368017112855051e-1 + 0.510878114959572e-1im
    t2 = 0.337315741065416 + 0.335449415919309im

    fr1 = (t1 - τ)*log(t1 - τ) + (t1 + τ)*log(t1 + τ) - 2*t1*log(t1) - τ*τ/t1
    fr2 = (t2 - τ)*log(t2 - τ) + (t2 + τ)*log(t2 + τ) - 2*t2*log(t2) - τ*τ/t2
    g_mass = ice_g0(model,Δπ) - s0*T + Tt*real(r1*fr1 + r2*fr2) #J/Kg
    return g_mass*molecular_weight(model,z)
end

function ice_g0(model::Ice2009,Δπ)
    g00 = -0.632_020_233_335_886e6
    g01 = 0.655_022_213_658_955
    g02 = -0.189_369_929_326_131e-7
    g03 = 0.339_746_123_271_053e-14
    g04 = -0.556_464_869_058_991e-21
    return evalpoly(Δπ,(g00,g01,g02,g03,g04))
end

function ice_r2(model::Ice2009,Δπ)
    r20 = - 0.725974574329220e2 - 0.781008427112870e2im
    r21 = -0.557107698030123e-4 + 0.464578634580806e-4im
    r22 = 0.234801409215913e-10 - 0.285651142904972e-10im
    return evalpoly(Δπ,(r20,r21,r22))
end

p_scale(model::Ice2009,z) = 101325.0
T_scale(model::Ice2009,z) = 273.16

molecular_weight(model::Ice2009,z) = 18.015268e-3*sum(z)

export Ice2009