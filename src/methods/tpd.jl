function mixture_fugacity(model,p,T,z;phase = :unknown, threaded = true, vol0 = nothing)
    V = p
    f = zeros(@f(Base.promote_eltype),length(z))
    return mixture_fugacity!(f,model,p,T,z;phase,threaded,vol0)
end


function mixture_fugacity!(f,model,p,T,z;phase = :unknown, threaded = true, vol0 = nothing)
    fugacity_coefficient!(f,model,p,T,z;phase,threaded,vol0)
    f .= f .* p .* z
    return f
end

function VT_mixture_fugacity!(f,model,V,T,z,p=pressure(model,V,T,z))
    f = VT_fugacity_coefficient!(f,model,V,T,z,p)
    f .= f .* p .* z
    return f
end

function VT_mixture_fugacity(model,V,T,z,p=pressure(model,V,T,z))
    f = zeros(@f(Base.promote_eltype),length(z))
    VT_mixture_fugacity!(f,model,V,T,z, p)
end

function tpd_obj(model, p, T, di, isliquid, cache = tpd_neq_cache(model,p,T,di,di), break_first = false)

    function f(α)
        w,vcache,dtpd,lnϕw,Hϕ = cache
        phase = isliquid ? :l : :v
        nc = length(model)
        w .= α .* α .* 0.25
        w ./= sum(w)
        volw0 = vcache[]
        lnϕw, volw = lnϕ!(lnϕw,model, p, T, w; phase=phase, vol0=volw0)
        dtpd .= log.(w) .+ lnϕw .- di
        fx = dot(w,dtpd) - sum(w) + 1
    end

    function g(df,α)
        w,vcache,dtpd,lnϕw,Hϕ = cache
        phase = isliquid ? :l : :v
        nc = length(model)
        w .= α .* α .* 0.25
        w ./= sum(w)
        lnϕw, volw = lnϕ!(lnϕw,model, p, T, w; phase=phase, vol0=vcache[])
        dtpd .= log.(w) .+ lnϕw .- di
        df .= dtpd .*  sqrt.(w)
        vcache[] = volw
        fx = dot(w,dtpd) - sum(w) + 1
        if fx < -1e-10 && break_first
            df .= 0
        end
        df
    end

    function fg(df,α)
        w,vcache,dtpd,lnϕw,Hϕ = cache
        phase = isliquid ? :l : :v
        nc = length(model)
        w .= α .* α .* 0.25
        w ./= sum(w)
        volw0 = vcache[]
        lnϕw, volw = lnϕ!(lnϕw,model, p, T, w; phase=phase, vol0=volw0)
        dtpd .= log.(w) .+ lnϕw .- di
        df .= dtpd .*  sqrt.(w)
        vcache[] = volw
        fx = dot(w,dtpd) - sum(w) + 1
        if fx < -1e-10 && break_first
            df .= 0
        end
        return fx,df
    end

    #=
    from thermopack:
    We see that ln Wi + lnφ(W) − di will be zero at the solution of the tangent plane minimisation.
    It can therefore be removed from the second derivative, without affecting the convergence properties.
    =#
    function fgh(df,d2f,α)
        w,vcache,dtpd,lnϕw,Hϕ = cache
        phase = isliquid ? :l : :v
        nc = length(model)
        w .= α .* α .* 0.25
        w ./= sum(w)
        volw0 = vcache[]
        lnϕw, ∂lnϕ∂nw, ∂lnϕ∂Pw, volw = ∂lnϕ∂n∂P(model, p, T, w, Hϕ; phase=phase, vol0=volw0)
        for i in 1:nc
            xi = sqrt(w[i])
            for j in 1:nc
                xj = sqrt(w[j])
                δij = Int(i == j)
                d2f[i,j] = δij + xi*xj*∂lnϕ∂nw[i,j]
            end
        end
        dtpd .= log.(w) .+ lnϕw .- di
        df .= dtpd .*  sqrt.(w)
        fx = dot(w,dtpd) - sum(w) + 1
        if fx < -1e-10 && break_first
            df .= 0
        end
        vcache[] = volw
        return fx,df,d2f
    end

    function h(d2f,α)
        w,vcache,dtpd,lnϕw,Hϕ = cache
        phase = isliquid ? :l : :v
        nc = length(model)
        w .= α.^2 ./ 4.0
        w ./= sum(w)
        volw0 = vcache[]
        lnϕw, ∂lnϕ∂nw, ∂lnϕ∂Pw, volw = ∂lnϕ∂n∂P(model, p, T, w, Hϕ; phase=phase, vol0=volw0)
        for i in 1:nc
            xi = sqrt(w[i])
            for j in 1:nc
                xj = sqrt(w[j])
                δij = Int(i == j)
                d2f[i,j] = δij + xi*xj*∂lnϕ∂nw[i,j]# + 0.5*αi*dtpd[i]
            end
        end
        vcache[] = volw
        d2f
    end

    obj = NLSolvers.ScalarObjective(f=f,g=g,fg=fg,fgh=fgh,h=h)
    optprob = OptimizationProblem(obj = obj,inplace = true)
end

function tpd_K0(model,p,T)
    return tp_flash_K0(model,p,T)
end

struct TPDKSolver end
struct TPDPureSolver end

"""
    wl,wv = tpd_solver(model,p,T,z,K0)

given p,T,z,K0, tries to perform a tangent-phase stability criterion with the given K value.
it tries both liquid and vapour phase. returns the resulting compositions `wl` and `wv`.
"""
function tpd_solver(model,p,T,z,K0,
    fz = mixture_fugacity(model,p,T,z),
    cache = tpd_cache(model,p,T,z),
    ss_solver = TPDKSolver(),
    isliquidz = true;
    break_first = false,
    lle = false,
    tol_trivial = 1e-5)

    vle = !lle
    ss_cache,newton_cache = cache
    _fz = ss_cache[2]::Vector
    _fz .= fz
    if any(isnan,fz)
        wl,wv = ss_cache[4],ss_cache[5]
        _0 = zero(eltype(fz))
        nan = _0/_0
        return wl,wv,nan,nan,nan,nan
    end

    #do initial sucessive substitutions
    stable_l, trivial_l, wl, vl = tpd_ss!(model,p,T,z,K0,true,ss_solver,ss_cache;tol_trivial)

    if isliquidz && vle #we suppose that vapour-vapour equilibria dont exist.
        stable_v, trivial_v, wv, vv = tpd_ss!(model,p,T,z,K0,false,ss_solver,ss_cache;tol_trivial)
    else
        stable_v, trivial_v, wv = true,true,ss_cache[5]
        vv = zero(eltype(wv))/zero(eltype(wv))
    end

    tpd_l = one(eltype(wl))*Inf
    tpd_v = one(eltype(wv))*Inf

    if trivial_l
        wl .= NaN
        tpd_l = first(wl)
    end

    if trivial_v
        wv .= NaN
        tpd_v = first(wv)
    end

    keep_going_l = !trivial_l && stable_l
    keep_going_v = !trivial_v && stable_v

    #fz = ϕ*p*z
    #log(fz) = log(p) + log(ϕi) + log(z)
    _fz .= log.(fz ./ p)
    #_fz .-= log.(z)
    di = _fz
    fxy = ss_cache[3]

    if !stable_l
        lnϕwl,vl = lnϕ!(fxy,model,p,T,wl,phase = :l)
        tpd_l = @sum(wl[i]*(lnϕwl[i] + log(wl[i]) - di[i])) - sum(wl) + 1
        tpd_l < 0 && break_first && return wl,wv,tpd_l,tpd_v,vl,vv
    end

    if !stable_v
        lnϕwv,vv = lnϕ!(fxy,model,p,T,wv,phase = :v)
        tpd_v = @sum(wv[i]*(lnϕwv[i] + log(wv[i]) - di[i])) - sum(wv) + 1
        tpd_v < 0 && break_first && return wl,wv,tpd_l,tpd_v,vl,vv
    end

    opt_options = OptimizationOptions(f_abstol = 1e-12,f_reltol = 1e-8)

    newton_cache[2][] = vl
    if keep_going_l
        α0l = 2 .* sqrt.(wl)
        prob_l = tpd_obj(model, p, T, di, true, newton_cache, break_first)
        res_l = Solvers.optimize(prob_l, α0l, LineSearch(Newton(linsolve = static_linsolve)), opt_options)
        αl = Solvers.x_sol(res_l)
        wl .= αl .* αl .* 0.25
        wl ./= sum(wl)
        tpd_l = Solvers.x_minimum(res_l)
        vl = newton_cache[2][]
        tpd_l < 0 && break_first && return wl,wv,tpd_l,tpd_v,vl,vv
    end


    #success,tpd_l_proposed = assert_correct_volume(fxy,model,p,T,wl,vl,:l,di)
    #!success && (tpd_l = tpd_l_proposed)

    newton_cache[2][] = vv
    if keep_going_v
        α0v = 2 .* sqrt.(wv)
        newton_cache[2][] = vv
        prob_v = tpd_obj(model, p, T, di, false, newton_cache, break_first)
        res_v = Solvers.optimize(prob_v, α0v, TrustRegion(Newton(linsolve = static_linsolve),Dogleg()), opt_options)
        αv = Solvers.x_sol(res_v)
        wv .= αv .* αv .* 0.25
        wv ./= sum(wv)
        vv = newton_cache[2][]
        tpd_v = Solvers.x_minimum(res_v)
    end

    #success,tpd_v_proposed = assert_correct_volume(fxy,model,p,T,wv,vv,:v,di)
    #!success && (tpd_v = tpd_v_proposed)
    return wl,wv,tpd_l,tpd_v,vl,vv
end

function tpd_cache(model,p,T,z)
    ss_cache = tpd_ss_cache(model,p,T,z)
    newton_cache = tpd_neq_cache(model,p,T,z)
    return ss_cache,newton_cache
end

function tpd_ss_cache(model,p,T,z)
    V = p
    K0 = similar(z,@f(Base.promote_eltype))
    fz,fxy = similar(K0),similar(K0)
    wl,wv = similar(K0),similar(K0)
    ss_cache = (K0,fz,fxy,wl,wv)
end

function tpd_neq_cache(model,p,T,z)
    V = p
    w = similar(z,@f(Base.promote_eltype))
    w,dtpd,lnϕw = similar(w),similar(w),similar(w)
    vcache = Base.RefValue{eltype(w)}(NaN)
    Hϕ = ∂lnϕ_cache(model, p, T, z, Val{false}())
    newton_cache = (w,vcache,dtpd,lnϕw,Hϕ)
end


function tpd_ss!(model,p,T,z,K0,is_liquid,solver = TPDKSolver(),cache = tpd_ss_cache(model,p,T,z,K0);tol_equil = 1e-10, tol_trivial = 1e-5, maxiter = 30)
    _tpd_ss!(model,p,T,z,K0,solver,is_liquid,cache,tol_equil,tol_trivial,maxiter)
end

#tpd K-value solver
#a port from MultiComponentFlash.jl
function _tpd_ss!(model,p,T,z,K0,solver::TPDKSolver,is_liquid,cache,tol_equil, tol_trivial,maxiter)

    #phase of the trial composition
    phase = is_liquid ? :l : :v

    #is this a trivial solution?
    trivial = false
    S = 1.0
    iter = 0
    done = false
    K,fz,fw,wl,wv = cache
    w = is_liquid ? wl : wv
    K .= K0
    v = zero(eltype(w))/zero(eltype(w))
    while !done
        iter += 1
        S = zero(eltype(w))
        @inbounds for c in eachindex(w)
            w_i = is_liquid ? z[c]/K[c] : z[c]*K[c]
            w[c] = w_i
            S += w_i
        end
        @. w /= S
        v = volume(model,p,T,w,phase = phase,vol0 = v)
        VT_mixture_fugacity!(fw,model,v,T,w)
        R_norm = zero(eltype(K))
        K_norm = zero(eltype(K))
        @inbounds for c in eachindex(K)
            sfw = S*fw[c]
            R = is_liquid ? sfw/fz[c] : fz[c]/sfw
            K[c] *= R
            R_norm += (R-1)^2
            K_norm += log(K[c])^2
        end
        # Two convergence criteria:
        # - Approaching trivial solution (K-values are all 1)
        # - Equilibrium for a small amount of the "other" phase,
        #   the single-phase conditions are not stable.
        trivial = K_norm < tol_trivial || !isfinite(R_norm)
        converged = R_norm < tol_equil

        # Termination of loop
        ok = trivial || converged
        done = ok || iter == maxiter
        if done && !ok
            trivial = true #we need to keep iterating
        end

    end
    stable = trivial || S <= 1 + tol_trivial
    return (stable, trivial, w, v)
end
#TPD pure solver.
#used in thermopack
#because we use K values, TPDKSolver is used instead
#but this is the legacy tpd ss solver.
function _tpd_ss!(model,p,T,z,w0,solver::TPDPureSolver,is_liquid,cache,tol_equil, tol_trivial,maxiter)
    phase = is_liquid ? :l : :v
    #is this a trivial solution?
    trivial = false
    stable = true
    iter = 0
    done = false
    di,fz,lnϕw,wl,wv = cache
    w = is_liquid ? wl : wv
    w .= w0
    di .= log.(fz ./ p)
    v = zero(eltype(w))/zero(eltype(w))
    S = zero(eltype(w))
    while !done
        iter += 1
        lnϕw, v = lnϕ!(lnϕw,model, p, T, w; phase=phase, vol0=v)
        S = zero(eltype(w))
        for i in eachindex(w)
            wi = exp(di[i]-lnϕw[i])
            w[i] = wi
            S += wi
        end
        @. w /= S
        R_norm = zero(eltype(w))
        K_norm = zero(eltype(w))
        tpd = one(eltype(w))
        for i in eachindex(w)
            wi,lnϕwi = w[i],lnϕw[i]
            sfw = S*exp(lnϕwi)*p*wi
            R = is_liquid ? sfw/fz[i] : fz[i]/sfw
            R_norm += (R-1)^2
            K_norm += log(wi/z[i])^2
            tpd += wi*(log(wi) + lnϕwi - di[i] - 1)
        end
        # Two convergence criteria:
        # - Approaching trivial solution (K-values are all 1)
        # - Equilibrium for a small amount of the "other" phase,
        #   the single-phase conditions are not stable.
        trivial = K_norm < tol_trivial || !isfinite(R_norm)
        converged = R_norm < tol_equil #|| tpd < 0

        # Termination of loop
        ok = trivial || converged
        done = ok || iter == maxiter
        if done && !ok
            trivial = true #we need to keep iterating
        end
    end
    stable = trivial || S <= 1 + tol_trivial
    return (stable, trivial, w, v)
end
"""
    tpd(model,p,T,z;break_first = false,lle = false,tol_trivial = 1e-5, di = nothing)

Calculates the Tangent plane distance function (`tpd`). It returns:

- a vector with trial phase compositions where `tpd < 0`
- a vector with the `tpd` values
- a vector with symbols indicating the phase of the input composition
- a vector with symbols indicating the phase of the trial composition

It iterates over each two-phase combination, starting from pure trial compositions, it does succesive substitution, then Gibbs optimization.

If the vectors are empty, then the procedure couldn't find a negative `tpd`. That is an indication that the phase is (almost) surely stable.

"""
function tpd(model,p,T,z,cache = tpd_cache(model,p,T,z);reduced = false,break_first = false,lle = false,tol_trivial = 1e-5,strategy = :default, di = nothing)
    check_arraysize(model,z)
    if !reduced
        model_reduced,idx_reduced = index_reduction(model,z)
        zr = z[idx_reduced]
    else
        model_reduced = model
        idx_reduced = z .== z
        zr = z
    end
    result = _tpd(model_reduced,p,T,zr,cache,break_first,lle,tol_trivial,strategy,di)
    values,comps,phase_z,phase_w = result
    idx_by_tpd = sortperm(values)
    for i in idx_by_tpd
        comps[i] = index_expansion(comps[i],idx_reduced)
    end
    return comps,values[idx_by_tpd], phase_z[idx_by_tpd], phase_w[idx_by_tpd]
    #do index expansion and sorting here
end

function _tpd(model,p,T,z,cache = tpd_cache(model,p,T,z),break_first = false,lle = false,tol_trivial = 1e-5,strategy = :default, di = nothing)
    #step 0: initialize values

    if strategy == :default || strategy == :wilson
        K = tpd_K0(model,p,T) #normally wilson
    else
        K = zeros(Base.promote_eltype(model,p,T,z),length(z))
    end
    cond = (model,p,T,z)
    fz,phasez,v = tpd_input_composition(model,p,T,z,di,lle)
    isliquidz = is_liquid(phasez)
    vle = !lle
    values = zeros(eltype(fz),0)
    comps = fill(fz,0)
    phase_z = Symbol[]
    phase_w = Symbol[]
    result = values,comps,phase_z,phase_w
    isnan(v) && return result

    #we asked for lle, but the input phase itself is a vapour.
    lle && !isliquidz && return result
    lle = lle || !isliquidz #if we have a vapour phase, dont calculate additional ones.
    if strategy == :default || strategy == :wilson
        #step 1: wilson initial values  (and the inverse)
        proposed = tpd_solver(model,p,T,z,K,fz,cache,TPDKSolver(),isliquidz;break_first,lle,tol_trivial)
        add_to_tpd!(result,cond,proposed,isliquidz,tol_trivial)
        length(values) >= 1 && break_first && return result

        K .= 1 ./ K
        lle = lle | any(is_vapour,phase_w)
        proposed = tpd_solver(model,p,T,z,K,fz,cache,TPDKSolver(),isliquidz;break_first,lle,tol_trivial)
        add_to_tpd!(result,cond,proposed,isliquidz,tol_trivial)
        length(values) >= 1 && break_first && return result
    end

    #step 3: trying pure components

    if strategy == :pure || strategy == :default
        ids = sortperm(z) #we start with the compounds in the lowest amount first.
        for i in ids
            z_pure!(K,i)
            lle = lle | any(is_vapour,phase_w)
            proposed = tpd_solver(model,p,T,z,K,fz,cache,TPDPureSolver(),isliquidz;break_first,lle,tol_trivial)
            add_to_tpd!(result,cond,proposed,isliquidz,tol_trivial)
            length(values) >= 1 && break_first && return result
        end
    end

    return result
end

function tpd_input_composition(model,p,T,z,di,lle)
    if di == nothing
        vl = volume(model,p,T,z,phase = :l)
        vv = volume(model,p,T,z,phase = :v)
        
        if lle
            v = vl
            if vl ≈ vv
                Π = pip(model,v,T,z) #identify phase with pip
                isliquidz = Π <= 1.0
            else
                isliquidz = true
                phasez = :liquid
            end
        else
            idx,v,_ = volume_label((model,model),p,T,z,(vl,vv))
            if vl ≈ vv
                Π = pip(model,v,T,z) #identify phase with pip
                isliquidz = Π <= 1.0
            elseif idx == 1 #v = vl
                isliquidz = true
            elseif idx == 2 #v = vv
                isliquidz = false
            else
                # this should not be reachable. the only case is when v = NaN, in that case
                #we catch that later.
                isliquidz = false
            end
        end
        phasez = isliquidz ? :liquid : :vapour
        fz = VT_mixture_fugacity(model,v,T,z,p)
        return fz,phasez,v
    else
        fz = p .* exp.(di)
        phasez = :unknown
        v = zero(eltype(fz))
        return fz,phasez,v
    end
end

function z_pure!(K,i)
    K .= 0
    K[i] = 1
    K
end

function z_norm(z,w)
    z_norm = zero(Base.promote_eltype(z,w))
    for i in 1:length(z)
        z_norm += log(w[i]/z[i])^2
    end
    return z_norm
end

function add_to_tpd!(result,cond,proposed,is_liquid,tol_trivial = 1e-5)
    phase_zk = is_liquid ? :liquid : :vapour
    values,comps,phase_z,phase_w = result
    wl,wv,tpd_l,tpd_v,vl,vv = proposed
    proposed_l = (wl,tpd_l,vl,:liquid)
    proposed_v = (wv,tpd_v,vv,:vapour)
    #=if !isnan(tpd_l)
        @show (phase_zk,:liquid)
    end=#
    added_l = _add_to_tpd!(result,cond,proposed_l,tol_trivial)
    if added_l
        values,comps,phase_z,phase_w = result
        push!(phase_z,phase_zk)
    end
    if is_liquid
        #=if !isnan(tpd_v)
            @show (phase_zk,:vapour)
        end =#
        added_v = _add_to_tpd!(result,cond,proposed_v,tol_trivial)
        if added_v
            push!(phase_z,phase_zk)
        end
    else
        added_v = false
    end
    return added_l || added_v
end

function _add_to_tpd!(result,cond,proposed,phase,tol_trivial = 1e-5)
    w,tpd,v,phase = proposed
    values,comps,phase_z,phase_w = result
    model,T,p,z = cond
    isnan(tpd) && return false
    any(isnan,w) && return false
    tpd >= 0 && return false
    maximum(w) < 0 && return false
    z_norm(z,w) < tol_trivial && return false
    min_tpd = minimum(values,init = Inf*one(eltype(values))) |> abs
    for i in 1:length(comps)
        dz = z_norm(comps[i],w)
        dz < tol_trivial && return false
        abs(tpd - values[i]) < min_tpd*tol_trivial && return false
    end
    push!(values,tpd)
    push!(comps,deepcopy(w))
    push!(phase_w,phase)
    return true

end #with those checks, we can be sure that the new tpd is a different composition.

function __z_test(z)
    nc = length(z)
    z_test = fill(1e-3*one(eltype(z)),(Int64(nc*(1+(nc-1)/2)),nc))
    for i in 1:nc
        z_test[i,i] = 1.
    end
    k = nc+1
    for i in 1:nc-1
        for j in i+1:nc
            z_test[k,i] = 0.5
            z_test[k,j] = 0.5
            k += 1
        end
    end
    z_test = z_test .* z'
    z_test = z_test ./ sum(z_test;dims=2)
end

function suggest_K(model,p,T,z,pure = split_model(model))
    lnϕz,v = lnϕ(model,p,T,z,threaded = false)
    K = similar(lnϕz)
    di = similar(lnϕz)
    di .= lnϕz .+ log.(z)
    for i in 1:length(z)
        vl = volume(pure[i],p,T,phase = :l)
        vv = volume(pure[i],p,T,phase = :v)
        lnϕv = VT_lnϕ_pure(pure[i],vv,T,p)
        lnϕl = VT_lnϕ_pure(pure[i],vl,T,p)
        tpd_v = lnϕv - di[i]
        tpd_l = lnϕl - di[i]
        if tpd_l < 0 && tpd_v < 0
            K[i] = exp(lnϕl[1])/exp(lnϕv[1])
        elseif tpd_l < 0 && tpd_v >= 0
            K[i] = exp(lnϕl[1])/exp(lnϕz[i])
        elseif tpd_l >= 0 && tpd_v < 0
            K[i] = exp(lnϕz[i])/exp(lnϕv[1])
        else
            sat_x = initial_points_bd_T(pure[i],T)
            psat = first(sat_x)
            K[i] = psat / p
        end
    end
    return K
end

function K0_lle_init_cache(model::EoSModel,p,T)
    pure = split_model(model)
    μ_pure = gibbs_free_energy.(pure, p, T)
    return μ_pure
end

function K0_lle_init_cache(model::ActivityModel,p,T)
    return zeros(length(model))
end

function K0_lle_act_coefficient(model::EoSModel,p,T,z,μ_pure)
    μ_pure = K0_lle_init_cache(model::EoSModel,p,T)
    μ_mixt = chemical_potential(model, p, T, z)
    R̄ = Rgas(model)
    return exp.((μ_mixt .- μ_pure) ./ R̄ ./ T) ./z
end

function K0_lle_act_coefficient(model::ActivityModel,p,T,z,μ_pure)
    return activity_coefficient(model,p,T,z)
end

function K0_lle_init(model::EoSModel, p, T, z)
    nc = length(model)
    z_test = __z_test(z)
    ntest = length(@view(z_test[1,:]))
    μ_pure = K0_lle_init_cache(model::EoSModel,p,T)
    γ1 = K0_lle_act_coefficient(model,p,T,@view(z_test[1,:]),μ_pure)
    γ = fill(γ1,1)
    for i in 2:ntest
        push!(γ,K0_lle_act_coefficient(model,p,T,@view(z_test[i,:]),μ_pure))
    end
    γ2 = γ
    _0 = zero(eltype(γ1))
    err = Inf*one(_0)
    i_min = 0
    j_min = 0
    ϕi = similar(γ1)
    for i in 1:ntest
        z_test_i = @view(z_test[i,:])
        γi = γ[i]
        for j in i+1:ntest
            z_test_j = @view(z_test[j,:])
            γj = γ[i]
            ϕ = _0
            ϕi_finite = 0
            for k in 1:length(model)
                ϕik = (z_test_i[k] - z[k])/(z_test_i[k] - z_test_j[k])
                ϕi[k] = ϕik
                if isfinite(ϕik)
                    ϕi_finite += 1
                    ϕ += ϕik
                end
            end
            ϕ = ϕ/ϕi_finite

            err_ij = _0
            for k in 1:length(model)
                zi,zj = z_test_i[k],z_test_j[k]
                logγᵢzᵢ = log(γi[k] * zi)
                logγⱼzⱼ = log(γj[k] * zj)
                err_ij += logγᵢzᵢ
                err_ij -= logγⱼzⱼ
                err_ij += logγᵢzᵢ*(ϕ*zi + (1 - ϕ)*zj - z[k])
            end
            #err_ij = sum(log.(γ[i].*z_test[i,:]) .-log.(γ[j].*z_test[j,:])+log.(γ[i].*z_test[i,:]).*(ϕ*z_test[i,:]+(1-ϕ)*z_test[j,:].-z))
            if err_ij < err
                err = err_ij
                i_min,j_min = i,j
            end
        end
    end

    K0 = γ[i_min]./ γ[j_min]

    #=
    #extra step, reduce magnitudes if we aren't in a single phase
    g0 = dot(z, K0) - 1.
    g1 = 1. - sum(zi/Ki for (zi,Ki) in zip(z,K0))

    if g0 < 0 && g1 <= 0 #we need to correct the value of g0
        g0i = z .* K0
        kz,idx = findmax(g0i)
        #the maximum K value is too great
    elseif g1 > 0 && g0 >= 0 #we need to correct the value of g1
        g1i = z ./ K0

    else #bail out here?

    end =#
    return K0
end

export tpd