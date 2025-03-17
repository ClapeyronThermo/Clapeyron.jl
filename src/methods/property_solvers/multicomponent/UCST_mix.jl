"""
    UCST_pressure(model::EoSModel,T;v0=x0_UCST_pressure(model,T))
    UCST_mix(model::EoSModel,T;v0=x0_UCST_pressure(model,T))

Calculates the Upper critical solution point of a mixture at a given Temperature.

returns:
- UCST Pressure [`Pa`]
- volume at UCST Point [`m³`]
- molar composition at UCST Point
"""
function UCST_pressure(model::EoSModel,T; v0=nothing)
    if v0 === nothing
        v0 = x0_UCST_pressure(model,T)
    end
    x0 = vcat(v0[1],v0[2][1:end-1])
    f! = (F,x) -> Obj_UCST_mix(model, F, x[2:end], exp10(x[1]), T)
    r  = Solvers.nlsolve(f!,x0,LineSearch(Newton2(x0)))
    sol = Solvers.x_sol(r)
    !all(<(r.options.f_abstol),r.info.best_residual) && (sol .= NaN)
    z_c = FractionVector(sol[2:end])
    V_c = exp10(sol[1])
    p_c = pressure(model, V_c, T, z_c)
    return (p_c, V_c, z_c)
end

"""
    x0_UCST_pressure(model::EoSModel,T)

Initial point for `UCST_pressure(model,T)`.

Returns a tuple, containing:
- Base 10 logarithm initial guess for liquid composition `[m³]`
- Initial guess for molar fractions at UCST Point (default: equimolar)
"""
function x0_UCST_pressure(model::EoSModel,T)
    x0 = Fractions.zeros(length(model))
    p = p_scale(model,x0)
    V  = x0_volume(model,p,T,x0,phase = :l)
    return (log10(V),x0)
end

function Obj_UCST_mix(model::EoSModel,F,z,V,T)
    z = FractionVector(z)
    L,detM = mixture_critical_constraint(model,V,T,z)
    F[1] = L
    F[2] = detM
    return F
end

"""
    UCST_temperature(model::EoSModel,p,x;v0 = nothing)

Calculates the Upper critical solution point of a mixture at a given pressure.


inputs:
- model: EoS model
- p: pressure [`Pa`]
- v0 (optional): an initial guess,consisting of a tuple of initial temperature, volume and composition.

returns:
- UCST Temperature [`K`]
- volume at UCST Point [`m³`]
- molar composition at UCST Point
"""
function UCST_temperature(model::EoSModel,p;v0 = nothing)
    if v0 === nothing
        T,logvol0,w0 = x0_UCST_temperature(model,p)
    else
        T,logvol0,w0 = v0
    end
    vol0 = exp10(logvol0)
    imax,_ = findmax(w0)
    ηw = η_from_v(model, vol0, T, w0)
    _,idx_max = findmax(w0)
    x0 = vcat(log(T),ηw,deleteat(w0,idx_max)) #select component with highest fraction as pivot
    f! = (F,x) -> Obj_UCST_temperature(model, F, x, p, idx_max)
    r  = Solvers.nlsolve(f!,x0,LineSearch(Newton2(x0)))
    sol = Solvers.x_sol(r)
    !all(<(r.options.f_abstol),r.info.best_residual) && (sol .= NaN)
    T_c = exp(sol[1])
    z_c = FractionVector(sol[3:end],idx_max)
    V_c = v_from_η(model,sol[2],T_c,z_c)
    p_c = pressure(model, V_c, T, z_c)
    return (T_c, V_c, z_c)
end

function Obj_UCST_temperature(model, F, x, p, i)
    η,T = x[2],exp(x[1])
    z = FractionVector(@view(x[3:end]),i)
    V = v_from_η(model,η,T,z)
    L,detM = mixture_critical_constraint(model,V,T,z)
    F[1] = L
    F[2] = detM
    F[3] = (pressure(model,V,T,z) - p)/p
    return F
end

function x0_UCST_temperature(model::EoSModel,p)
    x0 = Fractions.zeros(length(model))
    T = T_scale(model,x0)
    V  = x0_volume(model,p,T,x0,phase = :l)
    return (T,log10(V),x0)
end

const UCST_mix = UCST_pressure
