function obj_melting_pressure(model::CompositeModel,T,vs,vl,ps,μs)
    solid = solid_model(model)
    fluid = fluid_model(model)
    return μp_equality1_p(solid,fluid,vs,vl,T,ps,μs)
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

Calculates the melting pressure of a `CompositeModel` containing a solid and fluid phase EoS, at a specified temperature `T`.
You can pass a tuple of initial values for the volumes `(vs0,vl0)`.

returns:
- `pm` is melting Pressure `[Pa]`
- `vs` is melting solid volume at specified temperature `[m³]`
- `vl` is melting liquid volume at specified temperature `[m³]`
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
    if method.v0 == nothing
        v0 = x0_melting_pressure(model,T)
    else
        v0 = method.v0
    end
    vs0 = v0[1]
    vl0 = v0[2]
    _0 = zero(vs0*vl0*T*oneunit(eltype(model)))
    nan = _0/_0
    fail = (nan,nan,nan)
    valid_input = check_valid_2ph_input(vs0,vl0,true,T)
    if !valid_input
        return fail
    end
    fluid = fluid_model(model)
    solid = solid_model(model)
    ps,μs = equilibria_scale(fluid)
    result,converged = try_2ph_pure_pressure(solid,fluid,T,vs0,vl0,ps,μs,method)
    if converged
        return result
    else
        return fail
    end
end

function x0_melting_pressure(model::CompositeModel,T)
    solid = solid_model(model)
    liquid = fluid_model(model)
    z = SA[1.0]
    p = p_scale(liquid,z)
    vs00 = x0_volume(solid,p,T,z,phase = :s)
    vl00 = x0_volume(liquid,p,T,z,phase = :l)
    #=
    strategy:
    quadratic taylor expansion for Helmholtz energy
    isothermal compressibility approximation for pressure
   =#
    ps,μs = equilibria_scale(liquid)
    return solve_2ph_taylor(solid,liquid,T,vs00,vl00,ps,μs)
end

function Obj_Mel_Temp(model::EoSModel, F, T, V_s, V_l,p,p̄,T̄)
    z = SA[1.0]
    eos_solid(V) = eos(model.solid,V,T,z)
    eos_fluid(V) = eos(model.fluid,V,T,z)
    A_l,Av_l = Solvers.f∂f(eos_fluid,V_l)
    A_s,Av_s =Solvers.f∂f(eos_solid,V_s)
    g_l = muladd(-V_l,Av_l,A_l)
    g_s = muladd(-V_s,Av_s,A_s)
    F1 = -(Av_l+p)/p̄
    F2 = -(Av_s+p)/p̄
    F3 = (g_l-g_s)/(R̄*T̄)
    return SVector(F1,F2,F3)
end

struct ChemPotMeltingTemperature{V} <: ThermodynamicMethod
    T0::Union{Nothing,V}
    v0::V
    check_triple::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function ChemPotMeltingTemperature(;v0 = nothing,
                                    T0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10000)

    return ChemPotMeltingTemperature(v0,T0,check_triple,f_limit,atol,rtol,max_iters)
end

"""
    pm,vs,vl = melting_temperature(model::CompositeModel,p;v0=x0_melting_pressure(model,T))

Calculates the melting temperature of a `CompositeModel` containing a solid and fluid phase EoS, at a specified pressure `p`.
You can pass a tuple of initial values for the volumes `(vs0,vl0)`.

returns:
- Melting Temperature `[K]`
- melting solid volume at specified pressure `[m³]`
- melting liquid volume at specified pressure `[m³]`
"""
function melting_temperature(model::CompositeModel,p;kwargs...)
    method = init_preferred_method(melting_temperature,model,kwargs)
    return melting_temperature(model,p,method)
end
function init_preferred_method(method::typeof(melting_temperature),model::CompositeModel{<:EoSModel,<:EoSModel},kwargs)
    ChemPotMeltingTemperature(;kwargs...)
end

function melting_temperature(model::CompositeModel,p,method::ThermodynamicMethod)
    p = p*p/p
    return melting_temperature_impl(model,p,method)
end

function melting_temperature_impl(model::CompositeModel,p,method::ChemPotMeltingTemperature)
    solid = solid_model(model)
    fluid = fluid_model(model)
    T̄ = T_scale(fluid)
    p̄ = p_scale(fluid)
    if method.v0 == nothing
        v0 = x0_melting_temperature(model,p)
    else
        v0 = method.v0
    end
    _1 = 
    V0 = SVector(v0[1],log(v0[2]),log(v0[3]))
    f!(F,x) = Obj_Mel_Temp(model,F,x[1],exp(x[2]),exp(x[3]),p,p̄,T̄)
    results = Solvers.nlsolve(f!,V0,TrustRegion(Newton(),Dogleg()),NEqOptions(method))
    x = Solvers.x_sol(results)
    vs = exp(x[2])
    vl = exp(x[3])
    Tfus = x[1]
    converged = check_valid_eq2(solid_model(model),fluid_model(model),p,vs,vl,Tfus)
    if converged
    return Tfus, vs, vl
    else
        nan = zero(Tfus)/zero(Tfus)
        return nan,nan,nan
    end
end

function x0_melting_temperature(model::CompositeModel,p)
    Tt,pt,vs0,vl0,_ = triple_point(model)
    solid,fluid = solid_model(model),fluid_model(model)
    K0 = -dpdT_saturation(solid,fluid,vs0,vl0,Tt)*Tt*Tt/pt  
    #Clausius Clapeyron
    #log(P/Ptriple) = K0 * (1/T - 1/Ttriple)
    Tinv = log(p/pt)/K0 + 1/Tt
    T0 =  1/Tinv
    return T0,vs0,vl0
end

#=
init of pressure-based iterative methods
used by gibbs-based models
=#

function g_and_v(model,p,T,v;phase = :unknown)
    v = volume(model,p,T,SA[1.0],phase = phase,vol0 = v)
    g = VT_gibbs_free_energy(model,v,T,SA[1.0])
    return g,v
end

function g_and_v(model::GibbsBasedModel,p,T,v;phase = :unknown)
    return 𝕘∂𝕘dp(model,p,T,SA[1.0])
end


struct IsoGibbsMeltingPressure{V} <: ThermodynamicMethod
    p0::V
    check_triple::Bool
    f_limit::Float64
    atol::Float64
    rtol::Float64
    max_iters::Int
end

function IsoGibbsMeltingPressure(;p0 = nothing,
                                    check_triple = false,
                                    f_limit = 0.0,
                                    atol = 1e-8,
                                    rtol = 1e-12,
                                    max_iters = 10000)

    return IsoGibbsMeltingPressure(p0,check_triple,f_limit,atol,rtol,max_iters)
end

function init_preferred_method(method::typeof(melting_pressure),model::CompositeModel{<:EoSModel,<:GibbsBasedModel},kwargs)
    return IsoGibbsMeltingPressure(kwargs...)
end

function init_preferred_method(method::typeof(melting_pressure),model::CompositeModel{<:GibbsBasedModel,<:GibbsBasedModel},kwargs)
    return IsoGibbsMeltingPressure(kwargs...)
end

function init_preferred_method(method::typeof(melting_pressure),model::CompositeModel{<:GibbsBasedModel,<:EoSModel},kwargs)
    return IsoGibbsMeltingPressure(kwargs...)
end

function melting_pressure_impl(model::CompositeModel,T,method::IsoGibbsMeltingPressure)
    _1 = one(Base.promote_eltype(model,T))
    if method.p0 == nothing
        v0 = x0_melting_pressure(model,T)
        p0 = v0[3]*_1
    else
        p0 = method.p0*_1
    end
    
    solid = solid_model(model)
    fluid = fluid_model(model)
    nan = _1*NaN
    p,lnp,vl,vs = p0,log(p0),nan,nan
    max_iters = method.max_iters
    for i in 1:max_iters
        gl,vl = g_and_v(fluid,p,T,vl,phase = :liquid)
        gs,vs = g_and_v(solid,p,T,vs,phase = :solid)

        #f(lp) = gs(exp(lp)) - gl(exp(lp))
        #df(lp) = vs(exp(lp))*
        f = gs - gl
        df = p*(vs - vl)
        dlnp = -f/df
        !isfinite(dlnp) && break
        lnp = lnp + dlnp
        converged,_ = Solvers.convergence(lnp,lnp + dlnp,method.atol,method.rtol)
        p = exp(lnp)
        converged && return p,vs,vl
    end

    return nan,nan,nan
end
