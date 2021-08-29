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
