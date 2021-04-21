## Pure saturation conditions solver


## Mixture saturation solver
function bubble_pressure(model, T, x; v0 =nothing)
    TYPE = promote_type(eltype(T),eltype(x))
    lb_v = lb_volume(model,x)
    ts = T_scales(model,x)
    pmix = p_scale(model,x)
    #ps = p_scales(model,x)
    #Mollerup K
    #k0 =   (ps ./ p) .* exp.(5.42 .* (1.0 .- (ts ./ t)))
    k0 = 10.0
    if v0 === nothing
        y0    = k0 .*x./(1 .+x.*(k0 .- 1))
        y0    = y0 ./sum(y0)
        X     = x
        v0    = [log10(lb_v/0.45),
                 log10(lb_v/1e-4)]
        append!(v0,y0[1:end-1])
    end
    
    f! = (F,z) -> Obj_bubble_pressure(model, F, T, exp10(z[1]), exp10(z[2]), x,@view(z[3:end]),ts,pmix)
    #j! = (J,z) -> Jac_bubble_pressure(model, J, T, exp10(z[1]), exp10(z[2]), x[i,:], z[3:end])
    r  =Solvers.nlsolve(f!,v0)
    v_l = exp10(r.info.zero[1])
    v_v = exp10(r.info.zero[2])
    y = FractionVector(r.info.zero[3:end])
    P_sat = pressure(model,v_l,T,x)
    return (P_sat, v_l, v_v, y)
end

function Obj_bubble_pressure(model, F, T, v_l, v_v, x, y,ts,ps)
    #y = append!(y,1-sum(y))
    #y = [y...,1-sum(y)]
    y = FractionVector(y) #julia magic, check misc.jl
    μ_l = vt_chemical_potential(model,v_l,T,x)
    μ_v = vt_chemical_potential(model,v_v,T,y)
    p_l = pressure(model,v_l,T,x)
    p_v = pressure(model,v_v,T,y)
    
    for i in 1:length(x)
        F[i] = (μ_l[i]-μ_v[i])/(R̄*ts[i])
    end
    F[end] = (p_l-p_v)/ps
end
