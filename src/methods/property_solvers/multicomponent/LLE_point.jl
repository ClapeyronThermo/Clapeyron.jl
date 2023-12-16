#general initial point

function x0_lle_init(model::EoSModel, p, T, z, z0 = nothing)
    nc = length(model)
    if z0 == nothing
        z_test = initial_candidate_fractions(z)
    else
        z_test = near_candidate_fractions(z0,0.25)
    end
    ntest = length(z_test)
    γz = zeros(ntest,nc)
    for i in 1:ntest
        γz[i,:] = activity_coefficient(model,p,T,z_test[i]) .* z_test[i]
    end
    γz0 = activity_coefficient(model,p,T,z)  .* z
    err = ones(ntest)*Inf
    for i in 1:ntest
        #divided by norm to penalize points too close of the initial point
        err[i] = dnorm(γz0,γz[i])/dnorm(z_test[i],z)
    end
    (val, idx) = findmin(err)
    return z_test[idx]
end

## LLE pressure solver

function x0_LLE_pressure(model::EoSModel,T,x,p0 = nothing)
    pure = split_model(model)
    sat = saturation_pressure.(pure,T)
    vi = getindex.(sat,2)
    vlx = dot(vi,x)
    if p0 == nothing
        p0x = pressure(model,vlx,T,x)
    else
        p0x = p0
    end
    w = x0_lle_init(model,p0x,T,x)
    w2 = x0_lle_init(model,p0x,T,x,w)
    vlw = dot(vi,w2)
    prepend!(w2,log10.((vlx,vlw)))
    return w2[1:end-1]
end

"""
    LLE_pressure(model::EoSModel, T, x; v0 = x0_LLE_pressure(model,T,x))

calculates the Liquid-Liquid equilibrium pressure and properties at a given temperature.

Returns a tuple, containing:
- LLE Pressure `[Pa]`
- liquid volume of composition `x₁ = x` at LLE Point [`m³`]
- liquid volume of composition `x₂` at LLE Point  [`m³`]
- Liquid composition `x₂`
"""
function LLE_pressure(model::EoSModel, T, x; v0 =nothing)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        error("There is no LLE for a pure component")
    end
    x_r = x[idx_r]
    pmix = p_scale(model_r,x_r)
    if v0 === nothing
        v0 = x0_LLE_pressure(model_r,T,x_r)
    end

    f! = (F,z) -> Obj_bubble_pressure(model_r, F, T, exp10(z[1]), exp10(z[2]), x_r,z[3:end],pmix)
    options = NLSolvers.NEqOptions(maxiter = 1000) #this should converge in very few iters
    #putting the limit here allows to faster bail-out in case of unsucessful iteration
    r  =Solvers.nlsolve(f!,v0,LineSearch(Newton()),options)
    sol = Solvers.x_sol(r)
    v_l = exp10(sol[1])
    v_ll = exp10(sol[2])
    xx_r = FractionVector(sol[3:end])
    P_sat = pressure(model_r,v_l,T,x_r)
    xx = zeros(length(model))
    xx[idx_r] = xx_r
    return (P_sat, v_l, v_ll, xx)
end

"""
    LLE_temperature(model::EoSModel, p, x; T0 = x0_LLE_temperature(model,p,x))

calculates the Liquid-Liquid equilibrium temperature and properties at a given pressure.

Returns a tuple, containing:
- LLE Pressure `[Pa]`
- liquid volume of composition `x₁ = x` at LLE Point [`m³`]
- liquid volume of composition `x₂` at LLE Point  [`m³`]
- Liquid composition `x₂`
"""
function LLE_temperature(model::EoSModel,p,x;v0=nothing)
    model_r,idx_r = index_reduction(model,x)
    if length(model_r)==1
        error("There is no LLE for a pure component")
    end
    x_r = x[idx_r]
    pmix = p_scale(model_r,x_r)
    if v0 === nothing
        v0 = x0_LLE_temperature(model_r,p,x_r)
    end
    nc = length(model_r)
    f!(F,z) = Obj_bubble_temperature(model_r, F, p, z[1], exp10(z[2]), exp10(z[3]), x_r, z[4:nc+2],pmix)
    options = NLSolvers.NEqOptions(maxiter = 1000) #this should converge in very few iters
    #putting the limit here allows to faster bail-out in case of unsucessful iteration
    r  =Solvers.nlsolve(f!,v0[1:nc+2],LineSearch(Newton()),options)
    sol = Solvers.x_sol(r)
    T   = sol[1]
    v_l = exp10(sol[2])
    v_ll = exp10(sol[3])
    xx_r = FractionVector(sol[4:end])
    xx = zeros(length(model))
    xx[idx_r] = xx_r
    return T, v_l, v_ll, xx
end

function x0_LLE_temperature(model::EoSModel,p,x)
    #xx = Fractions.neg(x)
    pure = split_model(model)
    sat = saturation_temperature.(pure,p)
    T0 = 0.92*minimum(getindex.(sat,1)) #TODO: LLE points cannot be determined by pure data alone
    v0 = x0_LLE_pressure(model,T0,x,p)
    vi = 1-sum(v0[3:end])
    push!(v0,vi)
    prepend!(v0,T0)
    return v0

end
#=
function presents_LLE(model,p,T)
    pure = split_model(model)
    g_pure = gibbs_free_energy.(pure,p,T)

    function mixing_gibbs(x1)
        z = FractionVector(x1)
        log∑z = log(sum(z))
        g_mix = gibbs_free_energy(model,p,T,z)
        for i in 1:length(z)
            g_mix -= z[i]*(g_pure[i] + R̄*T*(log(z[i]) - log∑z))
        end
        return g_mix
    end

    f0(x1) = ForwardDiff.derivative(mixing_gibbs,x1)
    prob = Roots.ZeroProblem(f0,(0.9))
    res = Roots.solve(prob)
    return (res,mixing_gibbs(res))
end=#