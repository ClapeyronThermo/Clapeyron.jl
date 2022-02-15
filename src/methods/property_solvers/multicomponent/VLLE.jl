## VLLE Solver
function Obj_VLLE_pressure(model::EoSModel, F, T, v_l, v_ll, v_v, x, xx, y,ts,ps)
    x   = FractionVector(x)
    y   = FractionVector(y)
    xx  = FractionVector(xx)#julia magic, check misc.jl
    μ_l   = VT_chemical_potential(model,v_l,T,x)
    μ_ll  = VT_chemical_potential(model,v_ll,T,xx)
    μ_v   = VT_chemical_potential(model,v_v,T,y)
    p_l   = pressure(model,v_l,T,x)
    p_ll  = pressure(model,v_ll,T,xx)
    p_v   = pressure(model,v_v,T,y)
    n_c   = length(model)
    for i in 1:n_c
        F[i] = (μ_l[i]-μ_v[i])/(R̄*ts[i])
        F[i+n_c] = (μ_ll[i]-μ_v[i])/(R̄*ts[i])
    end
    F[end-1] = (p_l-p_v)/ps
    F[end]   = (p_ll-p_v)/ps
    return F
end
"""
    VLLE_pressure(model::EoSModel, T; v0 = x0_LLE_pressure(model,T))

calculates the Vapor-Liquid-Liquid equilibrium pressure and properties of a binary mixture at a given temperature.

Returns a tuple, containing:
- VLLE Pressure `[Pa]`
- Liquid volume of composition `x₁` at VLLE Point [`m³`]
- Liquid volume of composition `x₂` at VLLE Point  [`m³`]
- Vapour volume of composition `y` at VLLE Point  [`m³`]
- Liquid composition `x₁`
- Liquid composition `x₂`
- Liquid composition `y`
"""
function VLLE_pressure(model::EoSModel, T; v0 =nothing)
    if v0 === nothing
        v0 = x0_VLLE_pressure(model,T)
    end
    v0[4]
    ts = T_scales(model)
    pmix = p_scale(model,collect(FractionVector(v0[4])))
    n = length(model)
    nx = length(model) -1 
    len = 2*(n+2)
    x0 = vcat(v0...)
    Fcache = zeros(eltype(T),len)
    f! = (F,z) -> Obj_VLLE_pressure(model, F, T, exp10(z[1]), exp10(z[2]), exp10(z[3]), z[4:(nx+3)], z[(nx+4):(2nx+3)], z[(2nx+4):end],ts,pmix)
    r  = Solvers.nlsolve(f!,x0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_ll = exp10(sol[2])
    v_v = exp10(sol[3])
    x = FractionVector(sol[4:(nx+3)])
    xx = FractionVector(sol[(nx+4):(2nx+3)])
    y = FractionVector(sol[(2nx+4):end])
    P_sat = pressure(model,v_v,T,y)
    return (P_sat, v_l, v_ll, v_v, x, xx, y)
end

function x0_VLLE_pressure(model::EoSModel, T)
    pure = split_model(model)
    sat  = saturation_pressure.(pure,T)
    y0 = Fractions.zeros(length(model))
    x0    = [0.75,0.25] #if we change this, VLLE_pressure (and temperature) can be switched to more than two components.
    xx0 = Fractions.neg(x0)
    v_li = getindex.(sat,2)
    v_vi = last.(sat)
    v_l0  = dot(x0,v_li)
    v_ll0 = dot(xx0,v_li)
    v_v0  = dot(y0,v_vi)
    return (log10(v_l0),log10(v_ll0),log10(v_v0),x0[1:end-1],xx0[1:end-1],y0[1:end-1])
end
"""
    VLLE_temperature(model::EoSModel, p; T0 = x0_LLE_temperature(model,p))

calculates the Vapor-Liquid-Liquid equilibrium temperature and properties of a binary mixture at a given temperature.

Returns a tuple, containing:
- VLLE temperature `[K]`
- Liquid volume of composition `x₁` at VLLE Point [`m³`]
- Liquid volume of composition `x₂` at VLLE Point  [`m³`]
- Vapour volume of composition `y` at VLLE Point  [`m³`]
- Liquid composition `x₁`
- Liquid composition `x₂`
- Liquid composition `y`
"""
function VLLE_temperature(model,p;T0 = nothing)
    if T0 === nothing
        T0 = x0_VLLE_temperature(model,p)
    end
    f(z) = Obj_VLLE_temperature(model,z,p)
    fT = Roots.ZeroProblem(f,T0)
    T = Roots.solve(fT,Roots.Order0())
    P_sat, v_l, v_ll, v_v, x, xx, y = VLLE_pressure(model,T)
    return T, v_l, v_ll, v_v, x, xx, y
end

function x0_VLLE_temperature(model,p)
   return 1.7*sum(T_scales(model))/length(model)
end

function Obj_VLLE_temperature(model,T,p)
    p̃,v_l,v_ll,xx = VLLE_pressure(model,T)
    return p̃-p
end

#=
  0.002896 seconds (2.08 k allocations: 335.297 KiB)
(54504.07966631277, 
6.461304272627602e-5, 9.36874408395517e-5, 0.045108311104382584, 
[0.6842913136788397, 0.31570868632116034], 
[0.27359159385170784, 0.7264084061482922], 
[0.5774677072292492, 0.4225322927707508])
=#