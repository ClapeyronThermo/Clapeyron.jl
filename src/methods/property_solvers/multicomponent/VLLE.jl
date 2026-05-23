function __x0_VLLE_pressure(model,T,z = [0.5,0.5])
    check_arraysize(model,z)
    p_edge,_,_ = edge_pressure(model,T)
    K = tp_flash_K0(model,p_edge,T,z)
    y = rr_flash_vapor(K,z,zero((eltype(K))))
    x = rr_flash_liquid(K,z,one(eltype(K)))
    K_lle = K0_lle_init(model,p_edge,T,x)
    x1 = rr_flash_vapor(K,x,zero((eltype(K))))
    x2 = rr_flash_liquid(K,x,one(eltype(K)))
    vy = volume(model,p_edge,T,y,phase = :v)
    vx1 = volume(model,p_edge,T,x1,phase = :l)
    vx2 = volume(model,p_edge,T,x2,phase = :l)
    return FlashResult(SVector(x1,x2,y),SVector(0.,0.,0.),SVector(vx1,vx2,vy),FlashData(p_edge,T)) 
end

## VLLE Solver
function Obj_VLLE_pressure(model::EoSModel, F, T, η_l, η_ll, η_v, _x, _xx, _y)
    x   = FractionVector(_x)
    y   = FractionVector(_y)
    xx  = FractionVector(_xx)
    v_l = v_from_η(model,η_l,T,x)
    v_ll = v_from_η(model,η_ll,T,xx)
    v_v = v_from_η(model,η_v,T,y)
    w = (x,xx,y)
    v = (v_l,v_ll,v_v)
    F = μp_equality(model, F, Tspec(T), v, w)
    return F
end

"""
    VLLE_pressure(model::EoSModel, T; v0 = x0_LLE_pressure(model,T))

Calculates the Vapor-Liquid-Liquid equilibrium pressure and properties of a binary mixture at a given temperature `T`.

Returns a tuple, containing:
- VLLE Pressure `[Pa]`
- Liquid volume of composition `x₁` at VLLE Point `[m³]`
- Liquid volume of composition `x₂` at VLLE Point `[m³]`
- Vapour volume of composition `y` at VLLE Point `[m³]`
- Liquid composition `x₁`
- Liquid composition `x₂`
- Vapour composition `y`
"""
function VLLE_pressure(model::EoSModel, T; v0 =nothing)
    if v0 === nothing
        v00 = x0_VLLE_pressure(model,T)
        w0 = vcat(v00...)
    else
        w0 = copy(v0)
    end    
    nx = length(model) - 1
    idx_x = 4:(nx+3)
    idx_xx = (nx+4):(2nx+3)
    idx_y = (2nx+4):length(w0)
    w0[1] = η_from_v(model,exp10(w0[1]),T,FractionVector(w0[idx_x]))
    w0[2] = η_from_v(model,exp10(w0[2]),T,FractionVector(w0[idx_xx]))
    w0[3] = η_from_v(model,exp10(w0[3]),T,FractionVector(w0[idx_y]))
    f! = (F,z) -> @inbounds Obj_VLLE_pressure(model, F, T,
    z[1], z[2], z[3],
    z[idx_x], z[idx_xx], z[idx_y])

    r  = Solvers.nlsolve(f!,w0,LineSearch(Newton2(w0)))
    sol = Solvers.x_sol(r)
    !__check_convergence(r) && (sol .= NaN)
    x = FractionVector(sol[idx_x])
    xx = FractionVector(sol[idx_xx])
    y = FractionVector(sol[idx_y])
    v_l = v_from_η(model,sol[1],T,x)
    v_ll = v_from_η(model,sol[2],T,xx)
    v_v = v_from_η(model,sol[3],T,y)
    P_sat = pressure(model,v_v,T,y)
    return (P_sat, v_l, v_ll, v_v, x, xx, y)
end

function x0_VLLE_pressure(model::EoSModel, T, z = [0.5,0.5])
    #=
    pure = split_pure_model(model)
    sat  = saturation_pressure.(pure,T)
    y0 = Fractions.zeros(length(model))
    x0    = [0.75,0.25] #if we change this, VLLE_pressure (and temperature) can be switched to more than two components.
    xx0 = Fractions.neg(x0)
    v_li = getindex.(sat,2)
    v_vi = last.(sat)
    v_l0  = dot(x0,v_li)
    v_ll0 = dot(xx0,v_li)
    v_v0  = dot(y0,v_vi) =#
    vlle0 = __x0_VLLE_pressure(model,p,z)
    vl0,vl10,v_v0 = vlle0.volumes
    x0,xx0,y0 = vlle0.compositions
    #T0 = temperature(result)
    return (log10(v_l0),log10(v_ll0),log10(v_v0),x0[1:end-1],xx0[1:end-1],y0[1:end-1])
end

"""
    VLLE_temperature(model::EoSModel, p; T0 = x0_LLE_temperature(model,p))

Calculates the Vapor-Liquid-Liquid equilibrium temperature and properties of a binary mixture at a given pressure `p`.

Returns a tuple, containing:
- VLLE temperature `[K]`
- Liquid volume of composition `x₁` at VLLE Point `[m³]`
- Liquid volume of composition `x₂` at VLLE Point `[m³]`
- Vapour volume of composition `y` at VLLE Point `[m³]`
- Liquid composition `x₁`
- Liquid composition `x₂`
- Vapour composition `y`
"""
function VLLE_temperature(model::EoSModel,p;v0=nothing)
    if v0 === nothing
        v0 = x0_VLLE_temperature(model,p)
        w0 = vcat(v0...)
    else
        w0 = copy(v0)
    end
    nx = length(model) - 1
    idx_x = 5:(nx+4)
    idx_xx = (nx+5):(2nx+4)
    idx_y = (2nx+5):length(w0)
    T0 = w0[1]
    w0[2] = η_from_v(model,exp10(w0[2]),T0,FractionVector(w0[idx_x]))
    w0[3] = η_from_v(model,exp10(w0[3]),T0,FractionVector(w0[idx_xx]))
    w0[4] = η_from_v(model,exp10(w0[4]),T0,FractionVector(w0[idx_y]))
    f! = (F,z) -> @inbounds Obj_VLLE_temperature(model, F, p, z[1],
    z[2], z[3], z[4], z[idx_x], z[idx_xx], z[idx_y])
    options = NLSolvers.NEqOptions(maxiter = 1000)
    r  = Solvers.nlsolve(f!,w0,LineSearch(Newton2(w0)),options)
    sol = Solvers.x_sol(r)
    !__check_convergence(r) && (sol .= NaN)
    T  = sol[1]
    x = FractionVector(sol[idx_x])
    xx = FractionVector(sol[idx_xx])
    y = FractionVector(sol[idx_y])
    v_l = v_from_η(model,sol[2],T,x)
    v_ll = v_from_η(model,sol[3],T,xx)
    v_v = v_from_η(model,sol[4],T,y)
    return T, v_l, v_ll, v_v, x, xx, y
end

function __x0_VLLE_temperature(model,p,z = [0.5,0.5])
    check_arraysize(model,z)
    T_edge,_,_ = edge_temperature(model,p,z)
    K = tp_flash_K0(model,p,T_edge,z)
    y = rr_flash_vapor(K,z,zero((eltype(K))))
    x = rr_flash_liquid(K,z,one(eltype(K)))
    K_lle = K0_lle_init(model,p,T_edge,x)
    x1 = rr_flash_vapor(K,x,zero((eltype(K))))
    x2 = rr_flash_liquid(K,x,one(eltype(K)))
    vy = volume(model,p,T_edge,y,phase = :v)
    vx1 = volume(model,p,T_edge,x1,phase = :l)
    vx2 = volume(model,p,T_edge,x2,phase = :l)
    return FlashResult(SVector(x1,x2,y),SVector(0.,0.,0.),SVector(vx1,vx2,vy),FlashData(p,T_edge)) 
end

function x0_VLLE_temperature(model::EoSModel,p,z = [0.5,0.5])
    #=pure = split_pure_model(model)
    sat  = saturation_temperature.(pure,p)
    y0 = Fractions.zeros(length(model))
    T0 = 0.95*minimum(getindex.(sat,1))
    x0    = [0.75,0.25] #if we change this, VLLE_pressure (and temperature) can be switched to more than two components.
    xx0 = Fractions.neg(x0)
    v_li = getindex.(sat,2)
    v_vi = last.(sat)
    v_l0  = dot(x0,v_li)
    v_ll0 = dot(xx0,v_li)
    v_v0  = dot(y0,v_vi) =#
    vlle0 = __x0_LLE_temperature(model,p,z)
    vl0,vl10,v_v0 = vlle0.volumes
    x0,xx0,y0 = vlle0.compositions
    T0 = temperature(result)
    return (T0,log10(v_l0),log10(v_ll0),log10(v_v0),x0[1:end-1],xx0[1:end-1],y0[1:end-1])
end

function Obj_VLLE_temperature(model::EoSModel, F, p, T, ηl, ηll, ηv, _x, _xx, _y)
    x   = FractionVector(_x)
    y   = FractionVector(_y)
    xx  = FractionVector(_xx)
    v_l = v_from_η(model,ηl,T,x)
    v_ll = v_from_η(model,ηll,T,xx)
    v_v = v_from_η(model,ηv,T,y)
    w = (x,xx,y)
    v = (v_l,v_ll,v_v)
    F = μp_equality(model, F, Pspec(p,T), v, w)
    return F
end
