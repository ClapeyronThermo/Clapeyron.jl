function obj_melting_pressure(model::CompositeModel,F,T,vs,vl,p̄,T̄)
    z = SA[1.0]
    eos_solid(V) = eos(model.solid,V,T,z)
    eos_fluid(V) = eos(model.fluid,V,T,z)
    A_l,Av_l = Solvers.f∂f(eos_fluid,vl)
    A_s,Av_s =Solvers.f∂f(eos_solid,vs)
    μl = muladd(-vl,Av_l,A_l)
    μs = muladd(-vs,Av_s,A_s)
    ps = - Av_s
    pl = - Av_l
    F[1] = (μs - μl)/R̄/T̄
    F[2] = (ps - pl)/p̄
    return F
end

struct ChemPotMeltingPressure{V} <: ThermodynamicMethod
    v0::V
    check_triple::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function ChemPotMeltingPressure(;v0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10000)

    return ChemPotMeltingPressure(v0,check_triple,f_limit,atol,rtol,max_iters)
end

"""
    pm,vs,vl = melting_pressure(model::CompositeModel,T;v0=x0_melting_pressure(model,T))

Calculates the melting pressure of a `CompositeModel` containing a solid and fluid phase EoS, at a specified pressure.
You can pass a tuple of initial values for the volumes `(vs0,vl0)`.

returns:
- Melting Pressure [`Pa`]
- melting solid volume at specified temperature [`m³`]
- melting liquid volume at specified temperature [`m³`]
"""
function melting_pressure(model::CompositeModel,T;kwargs...)
    method = init_preferred_method(melting_pressure,model,kwargs)
    return melting_pressure(model,T,method)
end
function init_preferred_method(method::typeof(melting_pressure),model::CompositeModel{<:EoSModel,<:EoSModel},kwargs)
    ChemPotMeltingPressure(;kwargs...)
end

function melting_pressure(model::CompositeModel,T,method::ThermodynamicMethod)
    T = T*T/T
    return melting_pressure_impl(model,T,method)
end

function melting_pressure_impl(model::CompositeModel,T,method::ChemPotMeltingPressure)
    T̄ = T_scale(model.fluid)
    p̄ = p_scale(model.fluid)
    if method.v0 == nothing
        v0 = x0_melting_pressure(model,T)
    else
        v0 = method.v0
    end
    V0 = vec2(log(v0[1]),log(v0[2]),T)
    f!(F,x) = obj_melting_pressure(model,F,T,exp(x[1]),exp(x[2]),p̄,T̄)
    results = Solvers.nlsolve(f!,V0,LineSearch(Newton()))
    x = Solvers.x_sol(results)
    vs = exp(x[1])
    vl = exp(x[2])
    pfus = pressure(model.fluid, vl, T)
    return pfus, vs, vl
end

function x0_melting_pressure(model::CompositeModel,T)
    solid = solid_model(model)
    liquid = fluid_model(model)
    z = SA[1.0]
    vs00 = x0_volume_solid(solid,T,z)
    vl00 = x0_volume_liquid(liquid,T,z)
    fs(_V) = eos(solid,_V,T,z)
    fl(_V) = eos(liquid,_V,T,z)
    as,das,d2as = Solvers.f∂f∂2f(fs,vs00)
    al,dal,d2al = Solvers.f∂f∂2f(fs,vs00)
    ps = -das
    pl = -dal
    #=
    strategy:
    quadratic taylor expansion for helmholtz energy
    isothermal compressibility aproximation for pressure
    =#
    p_scale,μ_scale =  scale_sat_pure(liquid)
    As(_V) = as + das*(_V - vs00) + 0.5*d2as*(_V - vs00)^2
    Al(_V) = al + dal*(_V - vl00) + 0.5*d2al*(_V - vl00)^2
    Ps(_V) = log(_V/vs00)*(-_V*d2as) + ps
    Pl(_V) = log(_V/vl00)*(-_V*d2al) + pl
    function f0(F,x)
        vs = exp(x[1])
        vl = exp(x[2])
        asi,ali = As(vs),Al(vl)
        psi,pli = Ps(vs),Pl(vl)
        μsi,μli = asi + psi*vs, ali + pli*vl
        F[1] = (μsi - μli)*μ_scale
        F[2] = (psi - pli)*p_scale
    end
    x0 = vec2(log(vs00),log(vl00),T)
    res = Solvers.nlsolve(f0,x0,LineSearch(Newton()))
    v = Solvers.x_sol(res)
    vs0 = exp(v[1])
    vl0 = exp(v[2])
    return (vs0,vl0)
end
