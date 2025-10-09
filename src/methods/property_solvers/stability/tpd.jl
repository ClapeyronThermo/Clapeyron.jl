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
        #if fx < -1e-10 && break_first
        #    df .= 0
        #end
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
        #if fx < -1e-10 && break_first
        #    df .= 0
        #end
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
        #if fx < -1e-10 && break_first
        #    df .= 0
        #end
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

struct TPDKSolver end
struct TPDPureSolver end

function _tpd_and_v!(fxy,model,p,T,w,di,phase = :l)
    lnϕw,v = lnϕ!(fxy,model,p,T,w,phase = phase)
    tpd = @sum(w[i]*(lnϕw[i] + log(w[i]) - di[i])) - sum(w) + 1
    return tpd,v
end


"""
    wl,wv = tpd_solver(model,p,T,z,K0)

given p,T,z,K0, tries to perform a tangent-phase stability criterion with the given K value.
It tries both liquid and vapour phase. Returns the resulting compositions `wl` and `wv`.
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
        tpd_l,vl = _tpd_and_v!(fxy,model,p,T,wl,di,:l)
        tpd_l < 0 && break_first && return wl,wv,tpd_l,tpd_v,vl,vv
    end

    if !stable_v
        tpd_v,vv = _tpd_and_v!(fxy,model,p,T,wv,di,:v)
        tpd_v < 0 && break_first && return wl,wv,tpd_l,tpd_v,vl,vv
    end

    opt_options = OptimizationOptions(f_abstol = 1e-12,f_reltol = 1e-8,maxiter = 100)

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

function _tpd_fz_and_v!(solver::TPDKSolver,fxy,model,p,T,w,v0,liquid_overpressure = false,phase = :l)
    v = volume(model,p,T,w,phase = phase,vol0 = v0)
    if isnan(v) && liquid_overpressure && is_liquid(phase)
        #michelsen recomendation: when doing tpd, sometimes, the liquid cannot be created at the
        #specified conditions. try elevating the pressure at the first iter.
        v = volume(model,1.2p,T,w,phase = phase)
    end
    VT_mixture_fugacity!(fxy,model,v,T,w)
    return fxy,v,true
end

function _tpd_fz_and_v!(solver::TPDPureSolver,fxy,model,p,T,w,v0,liquid_overpressure = false,phase = :l)
    lnϕw,v = lnϕ!(fxy,model,p,T,w,phase = phase,vol0 = v0)
    if isnan(v) && liquid_overpressure && is_liquid(phase)
        #michelsen recomendation: when doing tpd, sometimes, the liquid cannot be created at the
        #specified conditions. try elevating the pressure at the first iter.
        lnϕw,v = lnϕ!(fxy,model,1.2p,T,w,phase = phase)
    end
    return fxy,v,true
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
    liquid_overpressure = false #if the liquid overpressure strategy was used
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
        _,v,liquid_overpressure = _tpd_fz_and_v!(solver,fw,model,p,T,w,v,liquid_overpressure,phase)
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
function _tpd_ss!(model,p,T,z,w0,solver::TPDPureSolver,_is_liquid,cache,tol_equil, tol_trivial,maxiter)
    phase = _is_liquid ? :l : :v
    #is this a trivial solution?
    trivial = false
    stable = true
    iter = 0
    done = false
    liquid_overpressure = false
    di,fz,lnϕw,wl,wv = cache
    w = _is_liquid ? wl : wv
    w .= w0
    di .= log.(fz ./ p)
    v = zero(eltype(w))/zero(eltype(w))
    S = zero(eltype(w))
    while !done
        iter += 1
        _,v,liquid_overpressure = _tpd_fz_and_v!(solver,lnϕw,model,p,T,w,v,liquid_overpressure,phase)
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
            R = _is_liquid ? sfw/fz[i] : fz[i]/sfw
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
function tpd(model,p,T,n,cache = tpd_cache(model,p,T,n);reduced = false,break_first = false,lle = false,tol_trivial = 1e-5,strategy = :default, di = nothing)
    z = n ./ sum(n)
    check_arraysize(model,z)
    if !reduced
        model_reduced,idx_reduced = index_reduction(model,z)
        zr = z[idx_reduced]
    else
        model_reduced = model
        idx_reduced = trues(length(model))
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
        K = tp_flash_K0(model,p,T,z) #normally wilson
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
                isliquidz = is_liquid(VT_identify_phase(model,v,T,z))
            else
                isliquidz = true
                phasez = :liquid
            end
        else
            idx,v,_ = volume_label((model,model),p,T,z,(vl,vv))
            if vl ≈ vv
                isliquidz = is_liquid(VT_identify_phase(model,v,T,z))
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

function suggest_K(model,p,T,z,pure = split_pure_model(model),volatiles = FillArrays.fill(true,length(model)),crit = FillArrays.fill(nothing,length(model)))
    lnϕz,v = lnϕ(model,p,T,z,threaded = false)
    K = similar(lnϕz)
    di = similar(lnϕz)
    log∑z = log(sum(z))
    di .= lnϕz .+ log.(z) .- log∑z
    for i in 1:length(z)
        if !volatiles[i]
            K[i] = 0
        end
        vl = volume(pure[i],p,T,phase = :l)
        vv = volume(pure[i],p,T,phase = :v)
        lnϕv = VT_lnϕ_pure(pure[i],vv,T,p)
        lnϕl = VT_lnϕ_pure(pure[i],vl,T,p)
        tpd_v = lnϕv - di[i]
        tpd_l = lnϕl - di[i]
        if vl ≈ vv
            if is_liquid(VT_identify_phase(pure[i],vv,T,SA[1.0])) || isnan(vv)
                ps,_,_ = saturation_pressure(pure[i],T,crit_retry = false)
                if isnan(ps)
                    tpd_l = Inf*abs(tpd_l)
                end
                tpd_v = Inf*abs(tpd_v)
            elseif is_vapour(VT_identify_phase(pure[i],vl,T,SA[1.0])) || isnan(vl)
                tpd_l = Inf*abs(tpd_l)
            end
        end
        if tpd_l < 0 && tpd_v < 0
            K[i] = exp(lnϕl[1])/exp(lnϕv[1])
        elseif tpd_l < 0 && tpd_v >= 0
            K[i] = exp(lnϕl[1])/exp(lnϕz[i])
        elseif tpd_l >= 0 && tpd_v < 0
            K[i] = exp(lnϕz[i])/exp(lnϕv[1])
        else #=tpd_l >= 0 && tpd_v >= 0=#
            #if tpd_l > tpd_v
                #K[i] = exp(lnϕv[1])
            #else
                #K[i] = exp(lnϕl[1])
            #end
            sat_x = extended_saturation_pressure(pure[i],T,crit[i])
            K[i] = first(sat_x)/p
        end
    end
    return K
end

function K0_lle_init(model::EoSModel, p, T, z)
    comps,tpds,_,_ = tpd(model,p,T,z,lle = true, strategy = :pure, break_first = true)
    if length(comps) == 1
        w = comps[1]
        β = one(eltype(w))
        for i in 1:length(z)
            β = min(β,z[i]/w[i])
        end
        β = 0.5*β
        w2 = (z .- β .*w)/(1 .- β)
        w2 ./= sum(w2)
        K = w ./ w2
    else
        K = ones(eltype(eltype(comps)),length(z))
    end
    return K
end


"""


"""
function double_tangency_points(model,p,T,z)
    idx = Vector{Int}[]
    n = length(model)
    for i in 1:n
        for j in i+1:n
            push!(idx,[i,j])
        end
    end
    binaries =  split_model(model,idx)
    V = p
    F = @f(Base.promote_eltype)
    comps = Vector{F}[]
    ijpairs = Tuple{Int,Int}[]
    cache = tpd_cache(binaries[1],p,T,[0.5,0.5])
    zij = zeros(F,2)
    phasez = Symbol[]
    phasew = Symbol[]
    for (ij,modelij) in pairs(binaries)
        ij_pair = idx[ij]
        ijt = (ij_pair[1],ij_pair[2])
        ij1,ij2 = ijt
        zij[1] = z[ij1]
        zij[2] = z[ij2]
        zij ./= sum(zij)
        
        compsij,_,phasezij,phasewij = tpd(modelij,p,T,zij,cache,strategy = :pure,break_first = true)
        if length(compsij) != 0
            append!(comps,compsij)
            for i in 1:length(compsij)
                push!(phasez,phasezij[1])
                push!(phasew,phasewij[1])
                push!(ijpairs,ijt)
            end
        end
    end
    return comps,ijpairs,phasez,phasew
end

export tpd