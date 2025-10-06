function obj_sublimation_pressure(model::CompositeModel,T,vs,vv,p_scale,μ_scale)
    solid = solid_model(model)
    fluid = fluid_model(model)
    return μp_equality1_p(solid,fluid,vs,vl,T,p_scale,μ_scale)
end

struct ChemPotSublimationPressure{V} <: ThermodynamicMethod
    v0::V
    check_triple::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function ChemPotSublimationPressure(;v0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 100)

    return ChemPotSublimationPressure(v0,check_triple,f_limit,atol,rtol,max_iters)
end

"""
    psub,vs,vv = sublimation_pressure(model::CompositeModel,T;v0=x0_sublimation_pressure(model,T))

Calculates the sublimation pressure `psub` of a `CompositeModel` containing a solid and fluid phase EoS, at a specified temperature `T`.
You can pass a tuple of initial values for the volumes `(vs0,vv0)`.

returns:
- Sublimation Pressure `[Pa]`
- Sublimation solid volume at specified temperature `[m³]`
- Sublimation vapour volume at specified temperature `[m³]`
"""
function sublimation_pressure(model::CompositeModel,T;kwargs...)
    method = init_preferred_method(sublimation_pressure,model,kwargs)
    return sublimation_pressure(model,T,method)
end

function init_preferred_method(method::typeof(sublimation_pressure),model::CompositeModel{<:EoSModel},kwargs)
    init_preferred_method(method,model.solid,kwargs)
end

function init_preferred_method(method::typeof(sublimation_pressure),model::EoSModel,kwargs)
    ChemPotSublimationPressure(;kwargs...)
end

function sublimation_pressure(model,T,method::ThermodynamicMethod)
    single_component_check(sublimation_pressure,model)
    T = T*T/T
    return sublimation_pressure_impl(model,T,method)
end

function sublimation_pressure_impl(model::CompositeModel,T,method::ChemPotSublimationPressure)
    if method.v0 == nothing
        v0 = x0_sublimation_pressure(model,T)
    else
        v0 = method.v0
    end
    vs0 = v0[1]
    vv0 = v0[2]
    _0 = zero(vs0*vv0*T*oneunit(eltype(model)))
    nan = _0/_0
    fail = (nan,nan,nan)
    valid_input = _is_positive((vs0,vv0,T))
    if !valid_input
        return fail
    end
    fluid = fluid_model(model)
    solid = solid_model(model)
    ps,μs = equilibria_scale(fluid)
    result,converged = try_2ph_pure_pressure(solid,fluid,T,vs0,vv0,ps,μs,method)
    if converged
        return result
    else
        return fail
    end
end

function x0_sublimation_pressure(model,T)
    #we can suppose we are in a low pressure regime, we treat the solid as a liquid,
    #and apply the zero pressure approximation.
    solid = solid_model(model)
    fluid = fluid_model(model)
    R̄ = Rgas(solid)
    z = SA[1.0]
    if solid isa GibbsBasedModel
        p = p_scale(solid)
        k1,k2 = calculate_gibbs_reference_state(model)
        g,dg,d2g = gibbs2_expansion(solid,p,T)
        vs_at_0 = dg - d2g*p
        g_ig = gibbs_energy(idealmodel(fluid),p,T)
        dg = g + k1 + k2*T - g_ig
        P0 = p*exp(dg/(R̄*T))
    else
        vs_at_0 = volume(solid,0.0,T,phase = :s)
        ares = a_res(solid, vs_at_0, T, z)
        lnϕ_s0 = ares - 1 + log(R̄*T/vs_at_0)
        P0 = exp(lnϕ_s0)
    end
    
    vv0 = Rgas(fluid)*T/P0
    vs0 = vs_at_0
    return vs0,vv0,P0
end

#=
lnphi_s = ln_phi0*V*(p - psub)


=#

function Obj_Sub_Temp(model::EoSModel, F, T, V_s, V_v, p, p̄, T̄)
    z = SA[1.0]
    eos_solid(V) = eos(model.solid,V,T,z)
    eos_fluid(V) = eos(model.fluid,V,T,z)
    A_v,Av_v = Solvers.f∂f(eos_fluid,V_v)
    A_s,Av_s =Solvers.f∂f(eos_solid,V_s)
    g_v = muladd(-V_v,Av_v,A_v)
    g_s = muladd(-V_s,Av_s,A_s)
    F1 = -(Av_v+p)/p̄
    F2 = -(Av_s+p)/p̄
    F3 = (g_v-g_s)/(R̄*T̄)
    return SVector(F1,F2,F3)
end

struct ChemPotSublimationTemperature{V} <: ThermodynamicMethod
    T0::Union{Nothing,V}
    v0::V
    check_triple::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function ChemPotSublimationTemperature(;v0 = nothing,
                                    T0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10000)

    return ChemPotSublimationTemperature(v0,T0,check_triple,f_limit,atol,rtol,max_iters)
end

"""
    pm,vs,vv = sublimation_temperature(model::CompositeModel,p;v0=x0_sublimation_pressure(model,T))

Calculates the sublimation temperature of a `CompositeModel` containing a solid and fluid phase EoS, at a specified pressure `p`.
You can pass a tuple of initial values for the volumes `(vs0,vv0)`.

returns:
- Sublimation Temperature `[K]`
- Sublimation solid volume at specified pressure `[m³]`
- Sublimation vapour volume at specified pressure `[m³]`
"""
function sublimation_temperature(model::CompositeModel,p;kwargs...)
    method = init_preferred_method(sublimation_temperature,model,kwargs)
    return sublimation_temperature(model,p,method)
end

function init_preferred_method(method::typeof(sublimation_temperature),model::CompositeModel{<:EoSModel},kwargs)
    init_preferred_method(method,model.solid,kwargs)
end

function init_preferred_method(method::typeof(sublimation_temperature),model::EoSModel,kwargs)
    ChemPotSublimationTemperature(;kwargs...)
end

function sublimation_temperature(model::CompositeModel,p,method::ThermodynamicMethod)
    p = p*p/p
    return sublimation_temperature_impl(model,p,method)
end

function sublimation_temperature_impl(model::CompositeModel,p,method::ChemPotSublimationTemperature)
    solid = solid_model(model)
    fluid = fluid_model(model)
    T̄ = T_scale(fluid)
    p̄ = p_scale(fluid)
    if method.v0 == nothing
        v0 = x0_sublimation_temperature(model,p)
    else
        v0 = method.v0
    end
    _1 = oneunit(p*1.0*one(eltype(solid))*one(eltype(fluid)))
    V0 = SVector(v0[1]*_1,log(v0[2])*_1,log(v0[3])*_1)
    f!(F,x) = Obj_Sub_Temp(model,F,x[1],exp(x[2]),exp(x[3]),p,p̄,T̄)
    results = Solvers.nlsolve(f!,V0,TrustRegion(Newton(),Dogleg()),NEqOptions(method))
    x = Solvers.x_sol(results)
    vs = exp(x[2])
    vv = exp(x[3])
    Tfus = x[1]
    converged = check_valid_eq2(solid,fluid,p,vs,vv,Tfus)
    if converged
    return Tfus, vs, vv
    else
        nan = zero(Tfus)/zero(Tfus)
        return nan,nan,nan
    end
end

function x0_sublimation_temperature(model::CompositeModel,p)
    Tt,pt,vs0,vl0,vv0 = triple_point(model)
    solid,fluid = solid_model(model),fluid_model(model)

    if solid isa GibbsBasedModel || fluid isa GibbsBasedModel
        K0 = -dpdT_saturation_gibbs(solid,fluid,pt,Tt,phase1 = :solid,phase2 = :gas)*Tt*Tt/pt
    else
        K0 = -dpdT_saturation(solid,fluid,vs0,vl0,Tt)*Tt*Tt/pt
    end

    #Clausius Clapeyron
    #log(P/Ptriple) = K0 * (1/T - 1/Ttriple)
    Tinv = log(p/pt)/K0 + 1/Tt
    T0 =  1/Tinv
    vs = volume(model,p,T0,phase = :s)
    vs = volume(model,p,T0,phase = :v)
    return T0,vs0,vv0
end

struct IsoGibbsSublimationPressure{V} <: ThermodynamicMethod
    p0::V
    check_triple::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function IsoGibbsSublimationPressure(;p0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10000)

    return IsoGibbsSublimationPressure(p0,check_triple,f_limit,atol,rtol,max_iters)
end

function init_preferred_method(method::typeof(sublimation_pressure),model::CompositeModel{<:GibbsBasedModel},kwargs)
    return IsoGibbsSublimationPressure(kwargs...)
end

function init_preferred_method(method::typeof(sublimation_pressure),model::GibbsBasedModel,kwargs)
    return IsoGibbsSublimationPressure(kwargs...)
end

function sublimation_pressure_impl(model::CompositeModel,T,method::IsoGibbsSublimationPressure)
    _1 = one(Base.promote_eltype(model,T))
    if method.p0 == nothing
        v0 = x0_sublimation_pressure(model,T)
        p0 = v0[3]*_1
    else
        p0 = method.p0*_1
    end

    solid = solid_model(model)
    fluid = fluid_model(model)
    nan = _1*NaN
    p,lnp,vl,vs = p0,log(p0),nan,nan
    valid_input = _is_positive((p,T))
    !valid_input && return (nan,nan,nan)
    max_iters = method.max_iters
    k1,k2 = calculate_gibbs_reference_state(model)
    for i in 1:max_iters
        gl,vl = g_and_v(fluid,p,T,vl,phase = :vapor)
        gs,vs = g_and_v(solid,p,T,vs,phase = :solid)

        #f(lp) = gs(exp(lp)) - gl(exp(lp))
        #df(lp) = vs(exp(lp))*
        f = gs - gl + k1 + k2*T
        df = p*(vs - vl)

        dlnp = -f/df
        !isfinite(dlnp) && break
        lnp = lnp + dlnp
        converged,_ = Solvers.convergence(lnp,lnp + dlnp,method.atol,method.rtol)
        p = exp(lnp)
        if converged 
            !_is_positive((p,vl,vs)) && break
            return p,vs,vl
        end
    end

    return nan,nan,nan
end

struct IsoGibbsSublimationTemperature{V} <: ThermodynamicMethod
    T0::V
    check_triple::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function IsoGibbsSublimationTemperature(;T0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10000)

    return IsoGibbsSublimationTemperature(T0,check_triple,f_limit,atol,rtol,max_iters)
end

function init_preferred_method(method::typeof(sublimation_temperature),model::CompositeModel{<:GibbsBasedModel},kwargs)
    return IsoGibbsSublimationTemperature(kwargs...)
end

function init_preferred_method(method::typeof(sublimation_temperature),model::GibbsBasedModel,kwargs)
    return IsoGibbsSublimationTemperature(kwargs...)
end

function sublimation_temperature_impl(model::CompositeModel,p,method::IsoGibbsSublimationTemperature)
    _1 = one(Base.promote_eltype(model,p))
    if method.T0 == nothing
        v0 = x0_sublimation_temperature(model,p)
        T0 = v0[1]*_1
    else
        T0 = method.T0*_1
    end

    solid = solid_model(model)
    fluid = fluid_model(model)
    nan = _1*NaN
    T,vl,vs = T0,nan,nan
    valid_input = _is_positive((p,T))
    !valid_input && return (nan,nan,nan)
    max_iters = method.max_iters
    k1,k2 = calculate_gibbs_reference_state(model)
    for i in 1:max_iters
        gl,sl,vl = g_and_sv(fluid,p,T,vl,phase = :vapor)
        gs,ss,vs = g_and_sv(solid,p,T,vs,phase = :solid)
        f = gs - gl + k1 + k2*T
        df = -(ss - sl - k2)
        dT = -f/df
        !isfinite(dT) && break
        T = T + dT
        converged,_ = Solvers.convergence(T,T + dT,method.atol,method.rtol)
        if converged 
            !_is_positive((T,vs,vl)) && return (nan,nan,nan)
            return T,vs,vl
        end
    end

    return nan,nan,nan
end
