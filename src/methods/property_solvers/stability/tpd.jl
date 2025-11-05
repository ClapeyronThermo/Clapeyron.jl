function tpd_cache(model,p,T,z,k0 = z)
    TT = Base.promote_eltype(model,p,T,z,k0)
    x1,x2,x3,x4 = similar(z,TT),similar(z,TT),similar(z,TT),similar(z,TT)
    vcache = Base.RefValue{TT}(NaN)
    Hϕ = ∂lnϕ_cache(model, p, T, z, Val{false}())
    return x1,x2,x3,x4,vcache,Hϕ
end

function tpd_obj(model, p, T, di, phase, cache = tpd_cache(model,p,T,di), break_first = false)

    function f(α)
        w,dtpd,_,_,vcache,Hϕ = cache
        nc = length(model)
        w .= α .* α .* 0.25
        w ./= sum(w)
        volw0 = vcache[]
        lnϕw, volw = lnϕ!(Hϕ, model, p, T, w; phase=phase, vol0=volw0)
        dtpd .= log.(w) .+ lnϕw .- di
        fx = dot(w,dtpd) - sum(w) + 1
    end

    function g(df,α)
        w,dtpd,_,_,vcache,Hϕ = cache
        nc = length(model)
        w .= α .* α .* 0.25
        w ./= sum(w)
        lnϕw, volw = lnϕ!(Hϕ, model, p, T, w; phase=phase, vol0=vcache[])
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
        w,dtpd,_,_,vcache,Hϕ = cache
        nc = length(model)
        w .= α .* α .* 0.25
        w ./= sum(w)
        volw0 = vcache[]
        lnϕw, volw = lnϕ!(Hϕ, model, p, T, w; phase=phase, vol0=volw0)
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
        w,dtpd,_,_,vcache,Hϕ = cache
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
        w,dtpd,_,_,vcache,Hϕ = cache
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

function _tpd_sum!(cache,model,p,T,w,di,v)
    lnϕw,_ = _tpd_fz_and_v!(cache,model,p,T,w,nothing,false,:unknown,v)
    tpd = @sum(w[i]*(lnϕw[i] + log(w[i]) - di[i])) - sum(w) + 1
    return tpd
end

function _tpd_fz_and_v!(cache,model,p,T,w,vol0,liquid_overpressure = false,phase = :liquid,_vol = nothing)
    if _vol === nothing
        vol = volume(model,p,T,w;phase,vol0)
    else
        vol = _vol*one(Base.promote_eltype(model,p,T,w))
    end
    
    fxy,v = lnϕ!(cache,model,p,T,w;vol)
    overpressure = false
    if isnan(v) && liquid_overpressure && is_liquid(phase)
        overpressure = true
        #michelsen recomendation: when doing tpd, sometimes, the liquid cannot be created at the
        #specified conditions. try elevating the pressure at the first iter.
        fxy,v = lnϕ!(cache,model,1.2p,T,w,phase = phase)
    end
    return fxy,v,overpressure
end

"""
    w,tpd,vw = tpd_solver(model,p,T,z,w0,dz = nothing;
                phasew = :liquid,
                break_first = true,
                tol_trivial = 1e-5,
                tol_equil = 1e-10,
                it_ss = 30)

Given `p`,`T`,`z`,`w0` and `dz = logϕ(z) .+ log(z)`, tries to perform a tangent-phase stability criterion with the given an initial input composition.
Returns a Tuple, containing:
 - molar composition at minimum tangent plane distance
 - Tangent plane distance at minimum
 - Vapour molar volume at minimum tangent plane distance `[m³·mol⁻¹]`
 - 
"""
function tpd_solver(model,p,T,z,w0,
    dz = nothing,
    cache = tpd_cache(model,p,T,z,w0);
    phasew = :unknown,
    break_first = false,
    tol_trivial = 1e-5,
    tol_equil = 1e-10,
    it_ss = 30,lle = false)

    w,_,_,dzz,vcache,Hϕ = cache

    if dz == nothing
        dz_temp,_,_ = tpd_input_composition(model,p,T,z,lle,cache)
        dzz .= dz_temp
    else
        dzz .= dz
    end

    if any(isnan,dzz)
        _0 = zero(eltype(dzz))
        nan = _0/_0
        return w,nan,nan
    end
    phase = phasew
    #do initial sucessive substitutions
    maxiter = it_ss
    w,tpd,vw,status = tpd_ss!(model,p,T,z,w0,cache;tol_trivial,tol_equil,maxiter,phase)
    stable,trivial = status
    if trivial
        w .= NaN
        tpd = NaN*tpd
    end

    keep_going = !trivial && stable

    if !stable
        tpd < 0 && break_first && return w,tpd,vw
    end

    
    vcache[] = vw
    
    if keep_going
        w,tpd,vw = tpd_optimization(model,p,T,z,w0,dzz,cache,phasew)
    end
    #success,tpd_v_proposed = assert_correct_volume(fxy,model,p,T,wv,vv,:vapour,di)
    #!success && (tpd_v = tpd_v_proposed)
    return w,tpd,vw
end

function tpd_ss!(model,p,T,z,w0,cache = tpd_cache(model,p,T,z,K0);phase = :liquid,tol_equil = 1e-10, tol_trivial = 1e-5, maxiter = 30)
    _tpd_ss!(model,p,T,z,w0,phase,cache,tol_equil,tol_trivial,maxiter)
end

function tpd_optimization(model,p,T,z,w0,di,cache = tpd_cache(model,p,T,z,K0),phasew = :liquid)
    w,_,_,dzz,vcache,Hϕ = cache
    α0 = 2 .* sqrt.(w0)
    prob = tpd_obj(model, p, T, dzz, phasew, cache)
    lb,ub = similar(α0),similar(α0)
    lb .= 0
    ub .= Inf
    opt_options = OptimizationOptions(f_abstol = 1e-12,f_reltol = 1e-10,maxiter = 100)
    res = Solvers.optimize(prob, α0, LineSearch(Solvers.Newton2(α0),Solvers.BoundedLineSearch(lb,ub)), opt_options)
    α = Solvers.x_sol(res)
    w .= α .* α .* 0.25
    w ./= sum(w)
    tpd = Solvers.x_minimum(res)
    vw = vcache[]
    return w,tpd,vw
end

#TPD pure solver.
function _tpd_ss!(model,p,T,z,w0,phase,cache,tol_equil,tol_trivial,maxiter)
    TT = Base.promote_eltype(model,p,T,z,w0)
    #is this a trivial solution?
    trivial = false
    stable = true
    iter = 0
    done = false
    liquid_overpressure = false
    w,_,fz,di,vcache,Hϕ = cache
    fz .= exp.(di) .* p
    w .= w0
    w ./= sum(w0)
    tpd,S,S_norm,v = TT(Inf),TT(Inf),TT(Inf),TT(NaN)
    while !done
        iter += 1
        lnϕw,v,liquid_overpressure = _tpd_fz_and_v!(Hϕ,model,p,T,w,v,liquid_overpressure,phase)
        S_old = S
        S = zero(TT)
        K_norm,tm,dtm = one(TT),zero(TT),zero(TT)
        for i in eachindex(w)
            wi = exp(di[i]-lnϕw[i])
            w[i] = wi
            S += wi
            
        end
        S_norm = abs(S_old - S)
        w ./=  S
        tpd = one(TT)
        for i in eachindex(w)
            wi,lnϕwi = w[i],lnϕw[i]
            K_norm += log(wi/z[i])^2
            tpd += wi*(log(wi) + lnϕwi - di[i] - 1)
        end
        #@show w,S_norm,q
        # Two convergence criteria:
        # - Approaching trivial solution (K-values are all 1)
        # - Equilibrium for a small amount of the "other" phase,
        #   the single-phase conditions are not stable.
        trivial = K_norm < tol_trivial #|| !isfinite(R_norm)
        converged = S_norm < tol_equil #|| tpd < 0

        # Termination of loop
        ok = trivial || converged
        done = ok || iter == maxiter
        done = done && iter > 1
        if done && !ok
            trivial = true #we need to keep iterating
        end
    end
    #@show S,iter,S_norm
    stable = trivial || S <= 1 + tol_trivial
    return (w,tpd,v,(stable,trivial))
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
function tpd(model,p,T,n,cache = tpd_cache(model,p,T,n);reduced = false,break_first = false,lle = false,tol_trivial = 1e-5,strategy = :default, di = nothing, verbose = false)
    z = n ./ sum(n)
    check_arraysize(model,z)

    if !reduced
        model_reduced,idx_reduced = index_reduction(model,n)
    else
        model_reduced,idx_reduced = model,fill(true,length(model))
    end

    zr = z[idx_reduced]
    eq = lle ? :lle : :vle
    model_reduced_cached = __tpflash_cache_model(model_reduced,p,T,z,eq)

    result = _tpd(model_reduced_cached,p,T,zr,cache,break_first,lle,tol_trivial,strategy,di,verbose)
    values,comps,phase_z,phase_w = result
    idx_by_tpd = sortperm(values)
    for i in idx_by_tpd
        comps[i] = index_expansion(comps[i],idx_reduced)
    end
    return comps,values[idx_by_tpd], phase_z[idx_by_tpd], phase_w[idx_by_tpd]
    #do index expansion and sorting here
end

function _tpd(model,p,T,z,cache = tpd_cache(model,p,T,z),break_first = false,lle = false,tol_trivial = 1e-5,strategy = :default, di = nothing,verbose = false)
    #step 0: initialize values
    
    cond = (model,p,T,z)

    if di != nothing
        WW = Base.promote_eltype(model,p,T,z,di)
        dz = similar(di,WW)
        dz .= di
        phasez = :unknown
        v = zero(WW)
    end

    dz,phasez,v = tpd_input_composition(model,p,T,z,lle,cache)
    TT = Base.promote_eltype(model,p,T,z,dz)
    K = similar(z,TT)
    K .= 0
    w_test = similar(K)
    isliquidz = is_liquid(phasez)
    vle = !lle
    values = zeros(TT,0)
    comps = fill(dz,0)
    phase_z = Symbol[]
    phase_w = Symbol[]
    result = values,comps,phase_z,phase_w
    isnan(v) && return result

    #we asked for lle, but the input phase itself is a vapour.
    lle && !isliquidz && return result
    lle_yet = lle & isliquidz #if we have a vapour phase, dont calculate additional ones.

    #plan what strategies we are gonna use:
    id_test = strategy == :default || strategy == :ideal_gas
    K_test = strategy == :default || strategy == :wilson || strategy == :K_values
    K_test = K_test && !lle_yet
    pure_test = strategy == :default || strategy == :pure

    #create a list of what test composition strategies we use
    tpd_strategies = tpd_plan(z,isliquidz,lle,id_test,K_test,pure_test)

    for strategy in tpd_strategies
        w,phasew,skip = tpd_test_composition!(strategy,cond,w_test,K,dz,verbose)
        skip && continue
        lle_yet && is_vapour(phasew) && continue
        proposed = tpd_solver(model,p,T,z,w,dz,cache;break_first,tol_trivial,phasew)
        added = add_to_tpd!(result,cond,proposed,phasez,phasew,tol_trivial)
        verbose && @info """
        $(tpd_print_strategy(strategy))
              Test composition:    $w_test
              Final composition:   $(proposed[1])
              Final tpd:           $(proposed[2])            
              Added to solution:   $added
        """
        added && break_first && return result
        is_vapour(phasew) && (lle_yet = lle_yet | any(is_vapour,phase_w))
    end

    return result
end

function tpd_plan(z,is_liquidz,lle,id_test,K_test,pure_test)
    plan = Tuple{Symbol,Symbol,Int}[]
    
    if is_liquidz && id_test && !lle
        push!(plan,(:ideal_gas,:vapour,0))
    end

    if K_test
        if is_liquidz
            lle || push!(plan,(:K,:vapour,1))
            push!(plan,(:K,:liquid,1))
            push!(plan,(:K,:liquid,-1))
            lle || push!(plan,(:K,:vapour,-1))
        else
            push!(plan,(:K,:liquid,1))
            push!(plan,(:K,:liquid,-1))
        end
    end

    if pure_test
        ids = sortperm(z)
        if is_liquidz
            for i in 1:length(z)
                push!(plan,(:pure,:liquid,ids[i]))
            end
            for i in 1:length(z)
                lle || push!(plan,(:pure,:vapour,ids[i]))
            end
        else
            for i in 1:length(z)
                push!(plan,(:pure,:liquid,ids[i]))
            end
        end
    end
    #display(plan)
    return plan

end

function tpd_test_composition!(strategy,conds,w_test,K,dz,verbose)
    plan,phase,ix = strategy
    model,p,T,z = conds
    skip = false
    skip_k = all(==(-1),K)

    is_k_plan = (plan == :K || plan == :invK) 
    if skip_k && is_k_plan
        skip = true
    end

    if plan == :ideal_gas
        w_test .= exp.(dz)
        w_test ./= sum(w_test)
    elseif plan == :pure
        z_pure!(w_test,ix)
    
    elseif is_k_plan && !skip_k
        if all(iszero,K)
            K .= tp_flash_K0(model,p,T,z) 
        end

        all(>(0),K) && verbose && @info "trying K-values:      $K"
        Kmin,Kmax = extrema(K)

        if Kmin > 1.0 || Kmax < 1.0
            verbose && @info "K-values unable to generate phase-split, skipping K-value tests."
            K .= -1
            skip = true
        end

        if !skip
            if ix == -1
                K .= 1 ./ K
            elseif ix != 1
                K .= K .^ ix
            end
            rr_flash_vapor!(w_test,K,z,false)
            w_test ./= sum(w_test)

            if ix == -1
                K .= 1 ./ K
            elseif ix != 1
                K .= K .^ -ix
            end
        end
    end
    return w_test,phase,skip
end

function tpd_print_strategy(strategy)
    plan,phase,ix = strategy
    if plan == :ideal_gas
        res = "Strategy: ideal gas, test phase: $phase"
    elseif plan == :K
        if ix == 1
            res = "Strategy: K-values, test phase: $phase"
        elseif ix == -1
            res = "Strategy: inverse K-values, test phase: $phase"
        else
            res = "Strategy: (K-values)^$(ix), test phase: $phase"
        end
        
    elseif plan == :pure
        res = "Strategy: pure initial point, test phase: $phase"
    else
        res = ""
    end
end

function tpd_input_composition(model,p,T,z,lle,cache = tpd_cache(model,p,T,z,di)) 
    TT = Base.promote_eltype(model,p,T,z)
    vl = volume(model,p,T,z,phase = :liquid)
    vv = volume(model,p,T,z,phase = :vapour)

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
    dz,_ = _tpd_fz_and_v!(last(cache),model,p,T,z,nothing,false,:unknown,v)
    logn = log(sum(z))
    dz .+= log.(z)
    dz .-= logn
    return copy(dz),phasez,v
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

function cosine_norm(z,w)
    znorm = sqrt(dot(z,z))
    wnorm = sqrt(dot(w,w))
    zwnorm = dot(z,w)
    return zwnorm/(wnorm*znorm)
end

function add_to_tpd!(result,cond,proposed,phasez,phasew,tol_trivial = 1e-5)
    w,tpd,v = proposed
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
        dzc = cosine_norm(comps[i],w)
        dz < tol_trivial && return false
        abs(tpd - values[i]) < min_tpd*tol_trivial && return false
    end
    push!(values,tpd)
    push!(comps,deepcopy(w))
    push!(phase_w,phasew)
    push!(phase_z,phasez)
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

function suggest_K(model,p,T,z,pure = split_pure_model(model),cache = nothing)
    lnϕz,v = lnϕ(model,p,T,z,cache,threaded = false)
    K = similar(lnϕz)
    di = similar(lnϕz)
    log∑z = log(sum(z))
    di .= lnϕz .+ log.(z) .- log∑z
    for i in 1:length(z)
        vl = volume(pure[i],p,T,phase = :liquid)
        vv = volume(pure[i],p,T,phase = :vapour)
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
            sat_x = extended_saturation_pressure(pure[i],T)
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