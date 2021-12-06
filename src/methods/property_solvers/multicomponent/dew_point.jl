## dew pressure solver
function x0_dew_pressure(model::EoSModel,T,y)
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
    _0 = zero(T+first(y))
    nan = _0/_0 
    sat_nan = (nan,nan,nan)
    replaceP = ifelse.(T_c .< T,true,false)

    eachx = eachcol(Diagonal(ones(eltype(y),length(y))))
#     Bi = second_virial_coefficient.(model,T,eachx)
    #using P_B(2B) as a sat aproximation
    #z = 1 + B/v
    #P_B = RT/v(1+B/v)
    #P_B(2B) = -RT/2B(1-B/2B)
    #P_B(2B) = -0.25*RT/B
    sat = [if !replaceP[i] saturation_pressure(pure[i],T) else sat_nan end for i in 1:length(pure)]
    
    P_sat = [tup[1] for tup in sat]
    V_l_sat = [tup[2] for tup in sat]
    V_v_sat = [tup[3] for tup in sat]
#     P_Bi = @. -0.25*R̄*T/Bi
    #=xP0 = yP
    #dot(x,P0) = P
    P = dot(x,P0)
    =#
    P⁻¹ = zero(T)
    V0_l = zero(T)
    V0_v = zero(T)
    Pi   = zero(y)
    for i in 1:length(y)
        if !replaceP[i]
            Pi[i] = P_sat[i][1]
            P⁻¹+=y[i]/Pi[i]
            V0_v += y[i]*V_v_sat[i]
        else 
            Pi[i] = pressure(pure[i],V_c[i],T)
            P⁻¹+=y[i]/Pi[i]
            V0_v += y[i]*V_c[i]*1.2
        end
    end
    #@show P_Bi
    #P = dot(x,P_Bi)
    P = 1/P⁻¹
    x = @. y*P/Pi
    xsum = 1/∑(x)
    x    = x.*xsum
    
    for i in 1:length(y)
        if !replaceP[i]
            V0_l += x[i]*V_l_sat[i]
        else
            V0_l += x[i]*V_c[i]
        end
    end
    
    prepend!(x,log10.([V0_l,V0_v]))
    return x
end

function dew_pressure(model::EoSModel, T, y; v0 =nothing)
    TYPE = promote_type(eltype(T),eltype(y))
#     lb_v = lb_volume(model,x)
    ts = T_scales(model,y)
    pmix = p_scale(model,y)
    if v0 === nothing
        v0 = x0_dew_pressure(model,T,y)
    end
    len = length(v0[1:end-1])
    #xcache = zeros(eltype(x0),len)
    Fcache = zeros(eltype(v0[1:end-1]),len)
    f!(F,z) = Obj_dew_pressure(model, F, T, exp10(z[1]), exp10(z[2]), z[3:end],y,ts,pmix)
    r  =Solvers.nlsolve(f!,v0[1:end-1],LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_v = exp10(sol[2])
    x = FractionVector(sol[3:end])
    P_sat = pressure(model,v_v,T,y)
    return (P_sat, v_l, v_v, x)
end

function Obj_dew_pressure(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    x   = FractionVector(x) #julia magic, check misc.jl
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

function dew_temperature(model,p,y)
    f(z) = Obj_dew_temperature(model,z,p,y)
    pure = split_model(model)
    sat = saturation_temperature.(pure,p)
    Ti   = zero(y)
    for i ∈ 1:length(y)
        if isnan(sat[i][1])
            Tc,pc,vc = crit_pure(pure[i])
            g(x) = p-pressure(pure[i],vc,x,[1.])
            Ti[i] = Roots.find_zero(g,(Tc))
        else
            Ti[i] = sat[i][1]
        end
    end
    T = Roots.find_zero(f,sum(Ti)/length(y))
    p,v_l,v_v,x = dew_pressure(model,T,y)
    return T,v_l,v_v,x
end

function Obj_dew_temperature(model,T,p,y)
    p̃,v_l,v_v,y = dew_pressure(model,T,y)
    return p̃-p
end

