## dew pressure solver
function x0_dew_pressure(model::EoSModel,T,y)
    pure = split_model(model)
    crit = crit_pure.(pure)
    
    T_c = [tup[1] for tup in crit]
    V_c = [tup[3] for tup in crit]
    _0 = zero(T+first(y))
    nan = _0/_0 
    sat_nan = (nan,nan,nan)
    replaceP = ifelse.(T_c .< T,true,false)
    sat = [if !replaceP[i] saturation_pressure(pure[i],T) else sat_nan end for i in 1:length(pure)]
    
    P_sat = [tup[1] for tup in sat]
    V_l_sat = [tup[2] for tup in sat]
    V_v_sat = [tup[3] for tup in sat]
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
    return μp_equality(model::EoSModel, F, T, v_l, v_v, FractionVector(x), y ,ts,ps)
end

function dew_temperature(model::EoSModel,p,y,T₀=nothing)
    TT = promote_type(typeof(p),eltype(y))
    nan = TT(NaN)
    if T₀ === nothing
        T₀::TT = x0_dew_temperature(model,p,y)
    end
    x0 = x0_dew_pressure(model,T₀,y)
    x = FractionVector(x0[3:end-1])
    v_l = exp10(x0[1])
    v_v = exp10(x0[2])
    cache = Ref{Tuple{TT,TT,TT,FractionVector{TT,Vector{TT}}}}((T₀,v_l,v_v,x))
    __f(z) = Obj_dew_temperature(model,z,p,y,cache)
    
    fT = Roots.ZeroProblem(__f,T₀)
    T::TT = Roots.solve(fT)
    return cache[]
end

function Obj_dew_temperature(model,T,p,y,cache)
    last_result = cache[]
    x0 = collect(last(last_result))
    prepend!(x0,(log10(last_result[2]),log10(last_result[3])))
    p̃,v_l,v_v,x= dew_pressure(model,T,y,v0 = x0)  
    cache[] = (T,v_l,v_v,x)
    return p̃-p
end

function x0_dew_temperature(model,p::T1,y::AbstractVector{T2}) where {T1,T2}
    TT = promote_type(T1,T2)
    pure = split_model(model)
    n = length(y)
    T̄ = zero(TT)
    for i ∈ 1:length(y)
        Tsat,_,_ = saturation_temperature(pure[i],p)
        if isnan(Tsat) && dew_temperature_T0i
            T̄ +=dew_temperature_T0i(pure[i],p)
        else
            T̄ += Tsat
        end
    end
    T̄ /= n
    return T̄::TT
end

function dew_temperature_T0i(model,p)
    Tc,_,vc = crit_pure(model)
    g(T) = p - pressure(model,vc,T)
    gi = Roots.ZeroProblem(g,Tc)
    Roots.solve(gi)
end