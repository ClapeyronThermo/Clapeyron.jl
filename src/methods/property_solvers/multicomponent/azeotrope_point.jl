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
    if v0 === nothing
        v0 = x0_azeotrope_pressure(model,T)
    end
    _,vl0,vv0,_ = __x0_bubble_pressure(model,T,v0)
    ηl = η_from_v(model,vl0,T,v0)
    ηv = η_from_v(model,vv0,T,v0)
    w0 = vcat(ηl,ηv,v0[1:end-1])
    f! = (F,z) -> Obj_az_pressure(model, F, T, z[1], z[2], z[3:end])
    r  =Solvers.nlsolve(f!,w0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    y = FractionVector(sol[3:end])
    v_l = v_from_η(model,sol[1],T,y)
    v_v = v_from_η(model,sol[2],T,y)
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

function Obj_az_pressure(model::EoSModel, F, T, ηl, ηv, x)
    xx = FractionVector(x)
    v_l = v_from_η(model,ηl,T,xx)
    v_v = v_from_η(model,ηv,T,xx)
    v = (v_l,v_v)
    w = (xx,xx)
    return μp_equality(model, F, Tspec(T), v, w)
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
    if v0 === nothing
        v0 = x0_azeotrope_temperature(model,p)
    end
    T0,vl0,vv0,_ = __x0_bubble_temperature(model,p,v0)
    ηl = η_from_v(model,vl0,T0,v0)
    ηv = η_from_v(model,vv0,T0,v0)
    w0 = vcat(T0,ηl,ηv,v0[1:end-1])
    f!(F,z) = Obj_azeotrope_temperature(model, F, p, z[1], z[2], z[3], z[4:end])
    r  = Solvers.nlsolve(f!,w0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    T = sol[1]
    x = FractionVector(sol[4:end])
    v_l = v_from_η(model,sol[2],T,x)
    v_v = v_from_η(model,sol[3],T,x)
    return T, v_l, v_v, x
end


function Obj_azeotrope_temperature(model::EoSModel, F, p, T, ηl, ηv, x)
    xx = FractionVector(x)
    v_l = v_from_η(model,ηl,T,xx)
    v_v = v_from_η(model,ηv,T,xx)
    v = (v_l,v_v)
    w = (xx,xx)
    return μp_equality(model, F, Pspec(p,T), v, w)
end

function x0_azeotrope_temperature(model,p)
    n = length(model)
    return Fractions.zeros(n)
end
