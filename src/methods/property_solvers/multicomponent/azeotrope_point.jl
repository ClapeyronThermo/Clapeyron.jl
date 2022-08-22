"""
    azeotrope_pressure(model::EoSModel, T; v0 = x0_azeotrope_pressure(model,T))

calculates the azeotrope pressure and properties at a given temperature.
Returns a tuple, containing:
- Azeotrope Pressure `[Pa]`
- liquid volume at Azeotrope Point [`m³`]
- vapour volume at Azeotrope Point [`m³`]
- Azeotrope composition
"""
function azeotrope_pressure(model::EoSModel, T; v0 = nothing)
    ts = T_scales(model)
    if v0 === nothing
        v0 = x0_azeotrope_pressure(model,T)
    end
    pmix = p_scale(model,v0)
    w0 = x0_bubble_pressure(model,T,v0)
    len = length(w0) - 1
    Fcache = zeros(eltype(v0),len)
    f! = (F,z) -> Obj_az_pressure(model, F, T, exp10(z[1]), exp10(z[2]), z[3:end],z[3:end],ts,pmix)
    r  =Solvers.nlsolve(f!,w0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    y = FractionVector(sol[3:end])
    P_sat = pressure(model,v_l,T,y)
    return (P_sat, v_l, v_v, y)
end

"""
    x0_azeotrope_pressure(model::EoSModel,T)

Initial point for `azeotrope_pressure(model,T)`.

Returns a vector, containing the initial guess azeotrope composition at a given temperature. defaults to equimolar
"""
function x0_azeotrope_pressure(model,T)
    n = length(model)
    return Fractions.zeros(n)
end

function Obj_az_pressure(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    return μp_equality(model::EoSModel, F, T, v_l, v_v, FractionVector(x), FractionVector(y),ts,ps)
end

"""
    azeotrope_temperature(model::EoSModel, T; v0 = x0_bubble_pressure(model,T,[0.5,0.5]))

Calculates the azeotrope temperature and properties at a given pressure.
Returns a tuple, containing:
- Azeotrope Temperature `[K]`
- liquid volume at Azeotrope Point [`m³`]
- vapour volume at Azeotrope Point [`m³`]
- Azeotrope composition
"""
function azeotrope_temperature(model::EoSModel,p;v0=nothing)
    ts = T_scales(model)
    if v0 === nothing
        v0 = x0_azeotrope_pressure(model,p)
    end
    pmix = p_scale(model,v0)
    w0 = x0_bubble_temperature(model,p,v0)
    len = length(w0) - 1
    Fcache = zeros(eltype(v0),len)
    w0[4:end] = v0
    f!(F,z) = Obj_azeotrope_temperature(model, F, p, z[1], exp10(z[2]), exp10(z[3]), z[4:end],z[4:end],ts,pmix)
    r  = Solvers.nlsolve(f!,w0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    T   = sol[1]
    v_l = exp10(sol[2])
    v_v = exp10(sol[3])
    x = FractionVector(sol[4:end])
    return T, v_l, v_v, x
end

function Obj_azeotrope_temperature(model::EoSModel, F, p, T, v_l, v_v, x, y,ts,ps)
    F = μp_equality(model::EoSModel, F, T, v_l, v_v, FractionVector(x), FractionVector(y),ts,ps)
    F[end] = (pressure(model,v_v,T,FractionVector(y)) - p)/ps
    return F
end
