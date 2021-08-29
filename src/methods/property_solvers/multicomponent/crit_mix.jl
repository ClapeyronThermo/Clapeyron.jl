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
    if v0 === nothing
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
