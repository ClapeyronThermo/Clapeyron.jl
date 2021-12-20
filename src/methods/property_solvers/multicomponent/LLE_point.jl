## LLE pressure solver
function x0_LLE_pressure(model::EoSModel,T,x)
    xx = 1 .-x
    pure = split_model(model)
    eachx = eachcol(Diagonal(ones(eltype(x),length(x))))
    sat = saturation_pressure.(pure,T)
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

function LLE_temperature(model,p,x;T0=nothing)
    f(z) = Obj_LLE_temperature(model,z,p,x)
    Ti   = T_scales(model,x).*1.5
    T = Roots.find_zero(f,sum(Ti)/length(Ti))
    if T0===nothing
        Ti   = T_scales(model,x).*1.5
        T = Roots.find_zero(f,sum(Ti)/length(Ti))
    else
        T = Roots.find_zero(f,T0)
    end
    p,v_l,v_ll,xx = LLE_pressure(model,T,x)
    return T,v_l,v_ll,xx
end

function Obj_LLE_temperature(model,T,p,x)
    p̃,v_l,v_ll,xx = LLE_pressure(model,T,x)
    return p̃-p
end

