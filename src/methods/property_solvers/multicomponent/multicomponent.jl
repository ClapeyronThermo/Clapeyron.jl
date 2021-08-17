## Bubble pressure solver
function x0_bubble_pressure(model::EoSModel,T,x)
    #TODO
    #on sufficiently large temps, 
    #the joule-thompson inversion occurs
    #making the virial coeff positive
    #on those cases, use an strategy that supposes pure gas on that side
    #Pbi = inf
    #xi = 0

    #check each T with T_scale, if treshold is over, replace Pi with inf
    pure = split_model(model)
    crit = crit_pure.(pure)
    
    T_c = [tup[1] for tup in crit]
    V_c = [tup[3] for tup in crit]
    _0 = zero(T+first(x))
    nan = _0/_0 
    sat_nan = (nan,nan,nan)
    replaceP = ifelse.(T_c .< T,true,false)

    eachx = eachcol(Diagonal(ones(eltype(x),length(x))))
#     Bi = second_virial_coefficient.(model,T,eachx)
    #using P_B(2B) as a sat aproximation
    #z = 1 + B/v
    #P_B = RT/v(1+B/v)
    #P_B(2B) = -RT/2B(1-B/2B)
    #P_B(2B) = -0.25*RT/B
    sat = [if !replaceP[i] sat_pure(pure[i],T) else sat_nan end for i in 1:length(pure)]
    
    P_sat = [tup[1] for tup in sat]
    V_l_sat = [tup[2] for tup in sat]
    V_v_sat = [tup[3] for tup in sat]
#     P_Bi = @. -0.25*R̄*T/Bi
    #=xP0 = yP
    #dot(x,P0) = P
    P = dot(x,P0)
    =#
    P = zero(T)
    V0_l = zero(T)
    V0_v = zero(T)
    Pi   = zero(x)
    for i in 1:length(x)
        if !replaceP[i]
            Pi[i] = P_sat[i][1]
            P+=x[i]*Pi[i]
            V0_l += x[i]*V_l_sat[i]
        else 
            Pi[i] = pressure(pure[i],V_c[i],T)
            P+=x[i]*Pi[i]
            V0_l += x[i]*V_c[i]
        end
    end
    #@show P_Bi
    #P = dot(x,P_Bi)
    y = @. x*Pi/P
    ysum = 1/∑(y)
    y    = y.*ysum
    
    for i in 1:length(x)
        if !replaceP[i]
            V0_v += y[i]*V_v_sat[i]
        else
            V0_v += y[i]*V_c[i]*1.2
        end
    end
    
    prepend!(y,log10.([V0_l,V0_v]))
    return y
end

function bubble_pressure(model::EoSModel, T, x; v0 =nothing)
    TYPE = promote_type(eltype(T),eltype(x))
#     lb_v = lb_volume(model,x)
    ts = T_scales(model,x)
    pmix = p_scale(model,x)
    if v0 === nothing
        v0 = x0_bubble_pressure(model,T,x)
    end
    len = length(v0[1:end-1])
    #xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(v0[1:end-1]),len)
    f! = (F,z) -> Obj_bubble_pressure(model, F, T, exp10(z[1]), exp10(z[2]), x,z[3:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    y = FractionVector(sol[3:end])
    P_sat = pressure(model,v_l,T,x)
    return (P_sat, v_l, v_v, y)
end

function Obj_bubble_pressure(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    y   = FractionVector(y) #julia magic, check misc.jl
    μ_l = VT_chemical_potential(model,v_l,T,x)
    μ_v = VT_chemical_potential(model,v_v,T,y)
    p_l = pressure(model,v_l,T,x)
    p_v = pressure(model,v_v,T,y)
    for i in 1:length(x)
        F[i] = (μ_l[i]-μ_v[i])/(R̄*ts[i])
    end
    F[end] = (p_l-p_v)/ps
    return F
end

    #j! = (J,z) -> Jac_bubble_pressure(model, J, T, exp10(z[1]), exp10(z[2]), x[i,:], z[3:end])
    #=
    _y0 = collect(FractionVector(v0[3:end]))
    @show _y0
    @show _pl0 = pressure(model,exp10(v0[1]),T,x)
    @show _pv0 = pressure(model,exp10(v0[2]),T,_y0)
    _P = 0.5*(_pl0+_pv0)
    for _ in 1:3
        _y0,_P = rr_bubble_pressure_refine(model,x,_y0,_P,T)
    end
    =#
# function bubble_pressure_rr(model, T, x; P = 40000)
#     sol0 = x0_bubble_pressure(model,T,x)
#     vl0 = exp10(sol0[1])
#     vv0 = exp10(sol0[2])




#     @show y0 = collect(FractionVector(sol0[3:end]))
#     @show vl0 = volume(model,P,T,x,phase=:l)
#     @show vv0 = volume(model,P,T,y0,phase=:v)
#     @show μ_l = vt_chemical_potential(model,vl0,T,x)
#     @show μ_v = vt_chemical_potential(model,vv0,T,y0)
#     y1 =  μ_l ./ μ_v .* x
#     @show y1 = y1 ./ sum(y1)
#     @show pl0 = pressure(model,vl0,T,x)
#     @show pv0 = pressure(model,vv0,T,y1)

#     @show P = (pl0 - pv0)/(log(pl0) - log(pv0))

#     @show vl = volume(model,P,T,x,phase=:l)
#     @show vv = volume(model,P,T,y1,phase=:v)
#     #=
#     μ_l = vt_chemical_potential(model,vl,T,x)
#     μ_v = vt_chemical_potential(model,vv,T,y)
#     K = log.(μ_v) ./ log.(μ_l)
#     y = K .* x
#     y = y./sum(y)
#     pl = pressure(model,vl,T,x)
#     pv = pressure(model,vv,T,y)
#     P = (pl+pv)/2
#     @show vl,vv,y
#     return y,P
#     =#
# end

## LLE pressure solver
function x0_LLE_pressure(model::EoSModel,T,x)
    xx = 1 .-x
    
    pure = split_model(model)

    eachx = eachcol(Diagonal(ones(eltype(x),length(x))))

    sat = sat_pure.(pure,T)
    
    P_sat = [tup[1] for tup in sat]
    V_l_sat = [tup[2] for tup in sat]
    
    V0_l = sum(x.*V_l_sat)
    V0_ll = sum(xx.*V_l_sat)
    
    prepend!(xx,log10.([V0_l,V0_ll]))
    return xx
end

function LLE_pressure(model::EoSModel, T, x; v0 =nothing)
    TYPE = promote_type(eltype(T),eltype(x))
#     lb_v = lb_volume(model,x)
    ts = T_scales(model,x)
    pmix = p_scale(model,x)
    if v0 === nothing
        v0 = x0_LLE_pressure(model,T,x)
    end
    len = length(v0[1:end-1])
    #xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(v0[1:end-1]),len)
    f! = (F,z) -> Obj_bubble_pressure(model, F, T, exp10(z[1]), exp10(z[2]), x,z[3:end],ts,pmix)
    r  =Solvers.nlsolve(f!,v0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_ll = exp10(sol[2])
    xx = FractionVector(sol[3:end])
    P_sat = pressure(model,v_l,T,x)
    return (P_sat, v_l, v_ll, xx)
end

## VLLE Solver
function Obj_three_phase(model::EoSModel, F, T, v_l, v_ll, v_v, x, xx, y,ts,ps)
    x   = FractionVector(x)
    y   = FractionVector(y)
    xx  = FractionVector(xx)#julia magic, check misc.jl
    μ_l   = VT_chemical_potential(model,v_l,T,x)
    μ_ll  = VT_chemical_potential(model,v_ll,T,xx)
    μ_v   = VT_chemical_potential(model,v_v,T,y)
    p_l   = pressure(model,v_l,T,x)
    p_ll  = pressure(model,v_ll,T,xx)
    p_v   = pressure(model,v_v,T,y)
    n_c   = length(model.components)
    for i in 1:n_c
        F[i] = (μ_l[i]-μ_v[i])/(R̄*ts[i])
        F[i+n_c] = (μ_ll[i]-μ_v[i])/(R̄*ts[i])
    end
    F[end-1] = (p_l-p_v)/ps
    F[end]   = (p_ll-p_v)/ps
    return F
end

function three_phase(model::EoSModel, T; v0 =nothing)
#     TYPE = promote_type(eltype(T),eltype(x))
#     lb_v = lb_volume(model,x)
    if v0 === nothing
        v0 = x0_three_phase(model,T)
    end
    ts = T_scales(model,[v0[4],1-v0[4]])
    pmix = p_scale(model,[v0[4],1-v0[4]])
    len = length(v0[1:end])
    #xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(v0[1:end]),len)
    f! = (F,z) -> Obj_three_phase(model, F, T, exp10(z[1]), exp10(z[2]), exp10(z[3]), z[4], z[5], z[6],ts,pmix)
    r  = Solvers.nlsolve(f!,v0[1:end],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_ll = exp10(sol[2])
    v_v = exp10(sol[3])
    x   = sol[4]
    xx = sol[5]
    y  = sol[6]
    x = FractionVector(x)
    xx = FractionVector(xx)
    y = FractionVector(y)
    P_sat = pressure(model,v_l,T,x)
    return (P_sat, v_l, v_ll, v_v, x, xx, y)
end

function x0_three_phase(model::EoSModel, T)
    pure = split_model(model)
    sat  = sat_pure.(pure,T)
    x0    = [0.75,0.25]
    xx0   = [0.25,0.75]
    y0    = [0.5,0.5]
    v_l0  = sat[1][2]*x0[1]+sat[2][2]*x0[2]
    v_ll0 = sat[1][2]*xx0[1]+sat[2][2]*xx0[2]
    v_v0  = sat[1][3]*y0[1]+sat[2][3]*y0[2]
    return [log10(v_l0),log10(v_ll0),log10(v_v0),x0[1],xx0[1],y0[1]]
end

function Obj_crit_mix(model::EoSModel,F,z,V,T)
    f(x) = eos(model,V,T,x)
    H(x) = ForwardDiff.hessian(f,x)/8.134/T
    L(x) = det(H(x))
    dL(x) = ForwardDiff.gradient(L,x)
    M(x) = [H(x)[1:end-1,:];transpose(dL(x))]
    F[1] = L(z)
    F[2] = det(M(z))
    return F
end
    
function crit_mix(model::EoSModel,z;v0=nothing)
    if v0 == nothing
        v0 = x0_crit_mix(model,z)
    end  
    f! = (F,x) -> Obj_crit_mix(model, F, z, exp10(x[1]), x[2])
    r  = Solvers.nlsolve(f!,v0[1:end],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    T_c = sol[2]
    V_c = exp10(sol[1])
    p_c = pressure(model, V_c, T_c, z)
    return (T_c, p_c, V_c)
end

function x0_crit_mix(model::EoSModel,z)
    pure = split_model(model)
    crit = crit_pure.(pure)
    n_c  = length(pure)
    V_c  = sum(z[i]*crit[i][3] for i in 1:n_c)
    T_c  = prod(crit[i][1]^z[i] for i ∈ 1:n_c)
    return [log10(V_c),T_c]
end

export bubble_pressure, LLE_pressure, three_phase, crit_mix
