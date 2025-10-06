struct IAPWS06 <: GibbsBasedModel
    components::Vector{String}
    s0::Float64
    Rgas::Float64
    references::Vector{String}
end

"""
    IAPWS06 <: GibbsBasedModel
    IAPWS06()

## Input parameters

None

## Description

IAPWS (International Association for the Properties of Water and Steam) Pure water Model for ice Ih, 2009 update.

## References

1. Feistel, R., & Wagner, W. (2006). A new equation of state for H2O ice Ih. Journal of Physical and Chemical Reference Data, 35(2), 1021–1047. [doi:10.1063/1.2183324](https://doi.org/10.1063/1.2183324)
2. IAPWS R10-06 (2009). Revised Release on the Equation of State 2006 for H2O Ice Ih
"""
IAPWS06

function IAPWS06(components = ["water"];s0 = -0.332733756492168e4,Rgas = 8.314472, idealmodel = nothing,verbose = false,reference_state = nothing)
    _components = format_components(components)
    references = default_references(IAPWS06)
    return IAPWS06(_components,s0,Rgas,references)
end

Rgas(model::IAPWS06) = model.Rgas
default_references(::Type{IAPWS06}) = ["IAPWS R10-06(2009)","10.1063/1.2183324"]
is_splittable(model::IAPWS06) = false

function eos_g(model::IAPWS06,p,T,z)
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

function ice_g0(model::IAPWS06,Δπ)
    g00 = -0.632_020_233_335_886e6
    g01 = 0.655_022_213_658_955
    g02 = -0.189_369_929_326_131e-7
    g03 = 0.339_746_123_271_053e-14
    g04 = -0.556_464_869_058_991e-21
    return evalpoly(Δπ,(g00,g01,g02,g03,g04))
end

function ice_r2(model::IAPWS06,Δπ)
    r20 = - 0.725974574329220e2 - 0.781008427112870e2im
    r21 = -0.557107698030123e-4 + 0.464578634580806e-4im
    r22 = 0.234801409215913e-10 - 0.285651142904972e-10im
    return evalpoly(Δπ,(r20,r21,r22))
end

p_scale(model::IAPWS06,z) = 101325.0
T_scale(model::IAPWS06,z) = 273.16

molecular_weight(model::IAPWS06,z) = 18.015268e-3*sum(z)

function x0_iapws06_fus(T)
    #ice ih melting pressure is only valid from 251.165 <= T <= 273.16
    Θ = T/273.16
    a1 = 0.119539337e7
    a2 = 0.808183159e5
    a3 = 0.333826860e4
    b1 = 0.300000e1
    b2 = 0.257500e2
    b3 = 0.103750e3
    π = 1 + a1*(1-Θ^b1) + a2*(1-Θ^b2) + a3*(1-Θ^b3)
    p = π*611.657
    return p
end

function x0_iapws06_sub(T)
    Θ = T/273.16
    a1 = -0.212144006e2
    a2 = 0.273203819e2
    a3 = -0.610598130e1
    b1 = 0.333333333e-2
    b2 = 0.120666667e1
    b3 = 0.170333333e1
    logπ = (a1*Θ^b1 + a2*Θ^b2 + a3*Θ^b3)/Θ
    p = exp(logπ)*611.657
    return p
end

function x0_pressure(model::IAPWS06,V,T,z)
    #sublimation pressure correlation
    _1 = one(Base.promote_eltype(model,T,V,z))
    pt = 611.57*_1
    Vst = 1.9652101514483842e-5

    if 251.165 <= T <= 273.16
        #it could be sublimation or melting correlation
        #check if the volume is lower or higher than the triple point      
        if V > Vst
            p = x0_iapws06_sub(T)
        else
            p = x0_iapws06_fus(T)
        end
    elseif T < 251.165
        #sublimation curve
        p = max(x0_iapws06_sub(T),pt)
        
    else #there is no ice ih in this range of temperature, return p_scale
        p = p_scale(model)
    end

    for i in 1:20
        if volume(model,p,T) <= V
            return p
        end
        p *= 2
    end
    return p
end

function x0_melting_pressure(model::CompositeModel{<:EoSModel,IAPWS06},T)
    solid = solid_model(model)
    liquid = fluid_model(model)
    z = SA[1.0]
    p  = x0_iapws06_fus(T)
    vs = volume(solid,p,T)
    vl = x0_volume(liquid,p,T,z,phase = :l)
    return vs,vl,p
end

function x0_melting_temperature(model::CompositeModel{<:EoSModel,IAPWS06},p)
    solid = solid_model(model)
    liquid = fluid_model(model)
    z = SA[1.0]
    f(T) = x0_iapws06_fus(T) - p
    prob = Roots.ZeroProblem(f,273.15)
    T = Roots.solve(prob)
    vs = volume(solid,p,T)
    vl = x0_volume(liquid,p,T,z,phase = :l)
    return T,vs,vl
end

function x0_sublimation_pressure(model::CompositeModel{<:EoSModel,IAPWS06},T)
    solid = solid_model(model)
    liquid = fluid_model(model)
    z = SA[1.0]
    p  = x0_iapws06_sub(T)
    vs = volume(solid,p,T)
    vv = Rgas(fluid)*T/p
    return vs,vv,p
end

function x0_sublimation_temperature(model::CompositeModel{<:EoSModel,IAPWS06},p)
    solid = solid_model(model)
    liquid = fluid_model(model)
    z = SA[1.0]
    f(T) = x0_iapws06_sub(T) - p
    prob = Roots.ZeroProblem(f,273.15)
    T = Roots.solve(prob)
    vs = volume(solid,p,T)
    vv = x0_volume(liquid,p,T,z,phase = :v)
    return T,vs,vv
end

function gibbsmodel_reference_state_consts(ice::IAPWS06,water::EmpiricHelmholtzModel)
    return :zero,0.0,0.0,0.0
end

function gibbsmodel_reference_state_consts(water::IAPWS06)
    return :dH,101325.0,273.15,6010.0
end

component_list(water::IAPWS06) = ["water"]

export IAPWS06
