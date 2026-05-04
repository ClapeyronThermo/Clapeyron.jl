struct MultiPhaseTPFlash{T} <: TPFlashMethod
    K0::Union{Vector{T},Nothing}
    n0::Union{Vector{Vector{T}},Nothing}
    K_tol::Float64
    ss_iters::Int
    nacc::Int
    second_order::Bool
    full_tpd::Bool
    max_phases::Int
    phase_iters::Int
    verbose::Bool
end

"""
    MultiPhaseTPFlash(;kwargs...)

Method to solve non-reactive multiphase (`np` phases), multicomponent (`nc` components) flash problem.

The flash algorithm uses successive stability tests to find new phases [1], and then tries to solve the system via Rachford-Rice and succesive substitution for `nc * np * ss_iters` iterations.

If the Rachford-Rice SS fails to converge, it proceeds to solve the system via Gibbs minimization in VT-space using lnK-β-ρ as variables [3].

The algorithm finishes when SS or the Gibbs minimization converges and all resulting phases are stable.

If the result of the phase equilibria is not stable, then it proceeds to add/remove phases again, for a maximum of `phase_iters` iterations.

### Keyword Arguments:

- `K0` (optional), initial guess for the constants K.
- `x0` (optional), initial guess for the composition of phase x.
- `y0` (optional), initial guess for the composition of phase y.
- `n0` (optional), initial guess for all compositions. It can be a matrix or a vector of vectors.
- `K_tol = sqrt(eps(Float64))`, tolerance to stop the calculation (`norm(lnK,1) < K_tol`).
- `ss_iters = 6`, number of Successive Substitution iterations to perform.
- `nacc = 5`, accelerate successive substitution method every nacc steps. Should be an integer bigger than 3. Set to 0 for no acceleration.
- `second_order = true`, whether to solve the Gibbs energy minimization using the analytical hessian or not. If set to `false`, the Gibbs minimization will be done using L-BFGS.
- `full_tpd` = false, whether to start with a simple K-split or using an intensive TPD search first.
- `max_phases = typemax(Int)`, the algorithm stops if there are more than `min(max_phases,nc)` phases.
- `phase_iters = 20`, the maximum number of solve-add/remove-phase iterations.

## References
1. Thermopack - Thermodynamic equilibrium algorithms reimplemented in a new framework. (2020, September 08). https://github.com/thermotools/thermopack. Retrieved May 4, 2024, from https://github.com/thermotools/thermopack/blob/main/docs/memo/flash/flash.pdf
2.  Okuno, R., Johns, R. T. T., & Sepehrnoori, K. (2010). A new algorithm for Rachford-Rice for multiphase compositional simulation. SPE Journal, 15(02), 313–325. [doi:10.2118/117752-pa](https://doi.org/10.2118/117752-pa)
3.  Adhithya, T. B., & Venkatarathnam, G. (2021). New pressure and density based methods for isothermal-isobaric flash calculations. Fluid Phase Equilibria, 537(112980), 112980. [doi:10.1016/j.fluid.2021.112980](https://doi.org/10.1016/j.fluid.2021.112980)
"""
function MultiPhaseTPFlash(;
    K0 = nothing,
    x0 = nothing,
    y0 = nothing,
    n0 = nothing,
    K_tol = sqrt(eps(Float64)),
    ss_iters = 6,
    nacc = 5,
    second_order = true,
    full_tpd = false,
    max_phases = typemax(Int),
    flash_result = nothing,
    phase_iters = 20,
    verbose = false) #TODO: find a better value for this

    if flash_result !== nothing
        comps = flash_result.compositions
        np = numphases(flash_result)
        np < 2 && incorrect_np_flash_error(MultiPhaseTPFlash,flash_result)
        ∑n = sum(flash_result.fractions)
        n00 = Vector{eltype(comps[1])}[]
        for i in 1:np
            push!(n00,collect(comps[i]))
        end
        return MultiPhaseTPFlash(;K0 = nothing,x0 = nothing, y0 = nothing,n0 = n00,K_tol,ss_iters,nacc,second_order,max_phases,phase_iters)
    end

    if K0 == x0 == y0 == n0 == nothing #nothing specified
    #is_lle(equilibrium)
        T = Nothing
        n00 = nothing
    else
        if !isnothing(K0) & isnothing(x0) & isnothing(y0) & isnothing(n0) #K0 specified
            T = eltype(K0)
            n00 = nothing
        elseif isnothing(K0) & !isnothing(x0) & !isnothing(y0) & isnothing(n0) #x0, y0 specified
            T = Base.promote_eltype(x0, y0, 1.0)
            n00 = Vector{T}[]
            push!(n00,convert(Vector{T},x0))
            push!(n00,convert(Vector{T},y0))
        elseif isnothing(K0) & isnothing(x0) & isnothing(y0) & !isnothing(n0)
            if n0 isa Matrix
                T = Base.promote_eltype(n0, 1.0)
                n00 = Vector{T}[]
                for i in 1:size(n0,1)
                    push!(n00,convert(Vector{T},@view(n0[1,:])))
                end
            elseif n0 isa AbstractVector && eltype(n0) <: AbstractVector && length(n0) > 0
                T = Base.promote_eltype(n0[1], 1.0)
                n00 = Vector{T}[]
                for i in 1:length(n0)
                    push!(n00,convert(Vector{T},n0[i]))
                end
            else
                throw(error("invalid specification of n0"))
            end
        else
            throw(error("invalid specification of initial points"))
        end
    end
    #check for nacc
    if nacc in (1,2,3) || nacc < 0
        throw(error("incorrect specification for nacc"))
    end

    ss_iters < 0 && throw(error("incorrect specification for ss_iters"))
    phase_iters < 1 && throw(error("incorrect specification for phase_iters"))

    return MultiPhaseTPFlash{T}(K0,n00,K_tol,ss_iters,nacc,second_order,full_tpd,max_phases,phase_iters,verbose)
end

is_vle(::MultiPhaseTPFlash) = false
is_lle(::MultiPhaseTPFlash) = false
is_unknown(::MultiPhaseTPFlash) = true

export MultiPhaseTPFlash

function index_reduction(m::MultiPhaseTPFlash,idx::AbstractVector)
    K0,n0,K_tol,ss_iters,nacc,second_order,full_tpd,max_phases,phase_iters,verbose = m.K0,m.n0,m.K_tol,m.ss_iters,m.nacc,m.second_order,m.full_tpd,m.max_phases,m.phase_iters,m.verbose
    K0 !== nothing && (K0 = K0[idx])
    if n0 !== nothing
        for i in 1:length(n0)
            n0i = n0[i]
            n0[i] = n0i[idx]
        end
    end
    # Delegate type resolution to keyword constructor to avoid eltype(K0) on nothing
    return MultiPhaseTPFlash(;K0,n0,K_tol,ss_iters,nacc,second_order,full_tpd,max_phases,phase_iters,verbose)
end

function tpd_cache end

function tp_flash_multi_cache(model,p,T,z)
    _tpd_cache = tpd_cache(model,p,T,z)
    Hϕ = _tpd_cache[end]
    F_cache = zeros(Base.promote_eltype(model,p,T,z),length(z))
    F_cache2 = similar(F_cache)
    F_cache3 = similar(F_cache)
    f1,f2,f3 = _tpd_cache[1],_tpd_cache[2],_tpd_cache[3]
    F3,F4,F5,ΔF1,ΔF2,xdem,Fdem = similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache)
    dem_cache = (F3,F4,F5,ΔF1,ΔF2,xdem,Fdem)
    comps_cache = fill(similar(f1),1)
    found_tpd = fill(similar(f1),0)
    found_tpd_lnphi = fill(similar(f1),0)
    found_tpd_volumes = similar(f1,0)
    v_cache = similar(f1)
    bi_cache = similar(f1)
    result_cache = comps_cache,bi_cache,v_cache
    phase_cache = (_tpd_cache,found_tpd,found_tpd_lnphi,found_tpd_volumes)
    ss_cache = result_cache,F_cache,F_cache2,F_cache3,f1,f2,f3,dem_cache,Hϕ
    return phase_cache,ss_cache
end

function resize_cache!(cache,np)
    phase_cache,ss_cache = cache
    nc = length(ss_cache[5])
    result_cache,F_cache,F_cache2,F_cache3,f1,f2,f3,dem_cache,Hϕ = ss_cache
    comps2,bi2,volumes2 = result_cache
    np_old = length(comps2)
    resize!(comps2,np)
    resize!(bi2,np)
    resize!(volumes2,np)
    for i in (np_old+1):np
        comps2[i] = similar(f1)
    end

    for fi in (F_cache,F_cache2,F_cache3)
        resize!(fi,nc*np)
    end

    for fi in dem_cache
        resize!(fi,nc*np)
    end
    return cache
end

function tp_flash_impl(model,p,T,z,method::MultiPhaseTPFlash)
    cached_model = __tpflash_cache_model(model,p,T,z,:unknown)
    return tp_flash_multi(cached_model,p,T,z,method)
end

function tp_flash_multi(model,p,T,nn,options = MultiPhaseTPFlash())
    z = nn ./ sum(nn)
    max_iter = options.phase_iters
    verbose = options.verbose
    TT = Base.promote_eltype(model,p,T,z)
    cache = tp_flash_multi_cache(model,p,T,z)
    phase_cache,ss_cache = cache
    idx_vapour = Ref(0)

    if options.n0 !== nothing
        #n phases provided. calculate volumes and fractions
        n0 = options.n0
        n_phases = length(n0)
        comps = Vector{TT}[]
        for i in 1:length(n0)
            push!(comps,convert(Vector{TT},n0[i]))
        end
        βi = initial_beta!(comps,z)
        volumes = similar(βi)
        for i in 1:length(comps)
            volumes[i] = volume(model,p,T,comps[i])
            if idx_vapour[] == 0
                if is_vapour(identify_phase(model,p,T,comps[i],vol = volumes[i]))
                    idx_vapour[] = i
                end
            end
        end
        δn_add = true

    elseif options.full_tpd #calculate full tpd.
        comps,tpds,_,phase_w = tpd(model,p,T,z,strategy = :default)
        min_tpd_verb = length(tpds) > 0 ? minimum(tpds) : zero(eltype(tpds))*NaN
        verbose && @info "[MPFLASH] full_tpd: candidates = $(length(comps)), min_tpd = $min_tpd_verb"
        if isempty(comps)
            # No unstable phase found by TPD: we can bail out with this.
            comps = [Vector{TT}(z)]
            vz = volume(model,p,T,z)
            volumes = [vz]
            βi = [one(eltype(z))]
            δn_add = false
            verbose && @info "[MPFLASH] full_tpd: fallback to single-phase"
        else
            βi = initial_beta!(comps,z)
            volumes = similar(βi)
            n_phases = length(comps)
            for i in 1:length(comps)
                phase_i = phase_w[i]
                volumes[i] = volume(model,p,T,comps[i],phase = phase_i)
                if idx_vapour[] == 0
                    if is_vapour(phase_i)
                        idx_vapour[] = i
                    end
                end
            end
            δn_add = true
        end
    else #split model manually
        res0 = pt_flash_x0(model,p,T,z,options)
        set_idx_vapour!(idx_vapour,model,res0)
        comps = [res0.compositions[1],res0.compositions[2]]
        volumes = [res0.volumes[1],res0.volumes[2]]
        βi = [res0.fractions[1],res0.fractions[2]]
        δn_add = true
    end

    result = FlashResult(comps, βi, volumes,FlashData(p,T))
    _result = (result,idx_vapour)
    δn_remove = false
    done = false
    converged = false
    iter = 0
    ss_converged = false
    neq_converged = false

    #step 1: tpd: we check if there is any unstable phase.

    #δn_add = _add_phases!(model,p,T,z,_result,cache,options)
    if !δn_add && length(comps) == 1
        g0 = _multiphase_gibbs(model,result)
        return FlashResult(comps, βi, volumes, FlashData(p,T,g0))
    end
    gmix = NaN*one(eltype(volumes))
    #step 2: main loop,iterate flashes until all phases are stable
    while !done
        iter += 1

        if δn_add  || δn_remove#step 2.1: δn != 0. we add or remove phases and find new candidate ones.
            n_phases = length(comps)
            #sucessive substitution iteration
            cache = resize_cache!(cache,n_phases)
            ss_converged = tp_flash_multi_ss!(model,p,T,z,_result,ss_cache,options)
            #verbose && @info "succesive substitution with $n_phases phases done. convergence: $ss_converged"
            #gibbs minimization
            if !ss_converged && minimum(βi) > 0 && all(isfinite,βi)
                neq_converged,gmix = tp_flash_multi_neq!(model,p,T,z,_result,ss_cache,options)
            end

            verbose && @info "[MPFLASH] flash with $n_phases phases done:"
            verbose && @info (repr("text/plain",result))
            #add/remove phases
            δn_remove = _remove_phases!(model,p,T,z,_result,cache,options)
            if !δn_remove
                δn_add = _add_phases!(model,p,T,z,_result,cache,options)
            else
                verbose && @info "[MPFLASH] removing phases:"
                verbose && @info (repr("text/plain",result))
                δn_add = false
            end
            if δn_add
                verbose && @info "[MPFLASH] adding phases:"
                verbose && @info (repr("text/plain",result))
            end
            converged = neq_converged || ss_converged
            no_new_phases = !δn_add && !δn_remove
            converged = converged && no_new_phases
            converged && ss_converged && (gmix = _multiphase_gibbs(model,result,idx_vapour[]))
        end
        done = iter > max_iter
        done = done || converged
        done = done || any(!isfinite,result.fractions) || minimum(result.fractions) < 0
        done = done || isone(numphases(result))
    end
    # Final cleanup: if duplicate phases remain (e.g., K≈1 leading to identical phases),
    # attempt a last consolidation to avoid returning duplicated phases.
    if numphases(result) > 1
        removed = _remove_phases!(model,p,T,z,_result,cache,options)
        verbose && @info("[MPFLASH] final cleanup: removed duplicate/degenerate phases")
    end
    if !isfinite(gmix)
        gmix = _multiphase_gibbs(model,result)
    end
    return FlashResult(result.compositions,result.fractions,result.volumes, FlashData(p,T,gmix))
end

function neq_converged(model,p,T,z,result)
    #TODO: add convergence criteria
    return true
end

function tp_flash_multi_ss!(model,p,T,z,_result,ss_cache,options)
    verbose = options.verbose
    result_cache,x0,x,xx,f1,f2,f3,dem_cache,Hϕ = ss_cache
    flash_result, idx_vapour = _result
    comps,βi,volumes = flash_result.compositions,flash_result.fractions,flash_result.volumes
    nc = length(z)
    np = length(βi)
    flash_result_dem = FlashResult(result_cache[1],result_cache[2],result_cache[3],FlashData(p,T))
    _result_dem = (flash_result_dem,idx_vapour)

    #we sort so the biggest phase fraction is at the end
    idx = sortperm(βi)
    comps .= comps[idx]
    βi .= βi[idx]
    volumes .= volumes[idx]
    if idx_vapour[] != 0
        new_idx_vapour = idx[idx_vapour[]]
        idx_vapour[] = new_idx_vapour
    end

    #set the initial values to 0
    x0 .= 0
    x .= 0
    xx .= 0

    #fill F with data
    xnp = comps[np]
    β = viewn(x0,nc,np)
    β = view(β,1:np)
    β .= βi
    for i in 1:(np-1)
        lnKi = viewn(x0,nc,i)
        xi = comps[i]
        # robust log ratio to avoid 0/0 → NaN
        for j in 1:nc
            xij,xnpj = xi[j],xnp[j]
            if xij == 0 && xnpj == 0
                lnKi[j] = zero(eltype(lnKi)) #log(1)
            else
                lnKi[j] = log(xij / xnpj)
            end
        end
    end

    #options
    ss_iters = options.ss_iters
    K_tol = options.K_tol
    nacc = options.nacc

    max_iters = ss_iters*np*nc
    itacc = 0
    converged = false
    F3,F4,F5,ΔF1,ΔF2,xdem,Fdem = dem_cache
    xdem .= 0
    for i in 1:max_iters
        itacc += 1
        fixpoint_multiphase!(x, x0, model, p, T, z, _result, ss_cache, verbose)

        #a phase needs to be removed.
        findfirst_duplicate_phases(flash_result,false) != (0,0) && break

        #store values in cache
        if minimum(βi) < 0 || any(!isfinite,βi)
            break
        end

        #acceleration
        if itacc == (nacc - 2)
            F3 .= x
        elseif itacc == (nacc - 1)
            F4 .= x
        elseif itacc == nacc
            itacc = 0
            F5 .= x
            # acceleration using DEM (1 eigenvalues)
            xdem = dem!(xdem, F5, F4, F3,(ΔF1,ΔF2))
            copyto!(flash_result_dem,flash_result)
            fixpoint_multiphase!(Fdem, xdem, model, p, T, z, _result_dem, ss_cache, verbose)
            gibbs = modified_gibbs(model,flash_result;vapour_phase_index=idx_vapour[])
            gibbs_dem = modified_gibbs(model,flash_result_dem;vapour_phase_index=idx_vapour[])
            if gibbs_dem < gibbs
                copyto!(flash_result,flash_result_dem)
                x .= Fdem
            end
        end
        converged = dnorm(x0,x,1) < K_tol
        converged && break
        x0 .= x
    end

    return converged
end

function fixpoint_multiphase!(F, x, model, p, T, z, _result, ss_cache, verbose = false)
    #given constant K, calculate β
    multiphase_RR_β!(F, x, z, _result , ss_cache, verbose)
    #given constant β, calculate K
    #fail early if we reach a non valid result on β
    result, idx_vapour = _result
    if minimum(result.fractions) < 0 || any(!isfinite,result.fractions)
        return F
    end
    multiphase_RR_lnK!(F, x, model, p, T, z, _result, ss_cache)
    return F
end

function RR_t!(t,x,β)
    nc = length(t)
    np = div(length(x),nc)
    t .= 1
    for l in 1:(np - 1)
        Kl = viewn(x,nc,l)
        βl = β[l]
        for i in 1:nc
            #K at phase i
            t[i] += βl*expm1(Kl[i])
        end
    end
    return t
end

function RR_β_ai!(a,F,i,nc,np)
    for j in 1:(np - 1)
        Kj = viewn(F,nc,j)
        a[j] = 1 - Kj[i]
    end
    return a
end

function RR_t_x!(x,K,t,z)
    for i in 1:length(x)
        Ki,ti,zi = K[i],t[i],z[i]
        if isinf(Ki)
            x[i] = zi
        else
            x[i] = Ki * zi / ti
        end
    end
end

function ls_restricted(φ::P,λ) where P
    _d = φ.d
    _x = φ.z
    λmax = λ
    #=
    x -  λ*d = 0
    λ = x/d
    λmax = minimum(xi/di for i in eachindex(x))
    =#
    for i in 1:length(_x)
        λmax = min(λmax,_x[i]/_d[i])
    end
    return λmax
end

#constant K, calculate β
function multiphase_RR_β!(F, x, z, _result, ss_cache, verbose)
    _,_,_,F_cache,f1,f2,f3,dem_cache,Hϕ = ss_cache
    result, idx_vapour = _result
    βi =result.fractions
    nc = length(z)
    np = length(βi)
    βF = viewn(F,nc,np)
    #I've added a verbose option here, which can be converted into function parameters later

    #transforms lnK to K
    multiphase_lnK_K!(F_cache,x,np,nc)
    if np == 2
        K = viewn(F_cache,nc,1)
        # Use log-K norm to detect trivial split (K ≈ 1 multiplicatively)
        Ktol = sqrt(eps(eltype(K)))
        if any(!isfinite, K) || __log_K_norm(K,false,false) <= Ktol
            βsol = clamp(βi[1], 0.0, 1.0)
            verbose && @info "[MPFLASH] RR_β: using previous β due to logK ≈ 0; βprev = $βsol"
        else
            βsol = rachfordrice(K,z)
        end
        βF[1] = βsol #two-phase solution does not require solving a neq problem
        βF[2] = 1 - βsol
        βi[1] = βsol
        βi[2] = 1 - βsol
        if verbose && !(isfinite(βsol) && 0 <= βsol <= 1)
            verbose && @info "[MPFLASH] RR_β: invalid βsol = $βsol, K= $K, z= $z"
        end
        return F
    end
    deactivated_phase = Ref(0)
    f = RR_obj(x, z, _result, ss_cache, deactivated_phase)
    β0 = copy(βi)
    resize!(β0,np - 1)

    lb_β = similar(β0)
    ub_β = similar(β0)
    lb_β .= 0
    ub_β .= Inf
    ls = Solvers.BoundedLineSearch(lb_β,ub_β)
    βsol = similar(β0)
    βresult = Solvers.optimize(f,β0,LineSearch(Newton(linsolve = static_linsolve),ls))
    i_deact = deactivated_phase[]
    βsol .= Solvers.x_sol(βresult)
    ∑β = sum(βsol)
    for i in 1:(np-1)
        βF[i] = βsol[i]
        βi[i] = βsol[i]
    end
    if i_deact != 0
        βF[i_deact] = -1
        βi[i_deact] = -1
    end
    #update last phase fraction
    βF[np] = 1 - ∑β
    βi[np] = 1 - ∑β
    return F
end

function RR_deactivated_phase!(deactivated_phase,β)
    i0 = deactivated_phase[]
    i0 != 0 && return i0

    ix = something(findfirst(<(8.881784197001252e-16),β),0)
    if ix != 0
        deactivated_phase[] = ix
    end
    return ix
end

function RR_deactivate_phase_gh!(ix,df,d2f)
    ix == 0 && return nothing
    if df !== nothing
        df[ix] = 0
    end
    if d2f !== nothing
        j = size(d2f,1)
        d2f1 = @view d2f[ix,:]
        d2f2 = @view d2f[:,ix]
        d2f1 .= 0
        d2f2 .= 0
        d2f[ix,ix] = 1
    end
    return nothing
end

function RR_obj(x, z, _result, ss_cache, deactivated_phase)
    #okuno et al. (2010): objetive function, gradient, hessian
    function f(β)
        _,_,_,F_cache,f1,f2,f3,dem_cache,Hϕ = ss_cache
        result, idx_vapour = _result
        nc,np = length(z),numphases(result)
        t = RR_t!(f3,x,β) #we use lnK here.
        fx = zero(eltype(β))
        for i in 1:nc
            fx -= z[i]*log(abs(t[i]))
        end
        return fx
    end

    function g(df,β)
        _,_,_,F_cache,f1,f2,f3,dem_cache,Hϕ = ss_cache
        result, idx_vapour = _result
        nc,np = length(z),numphases(result)
        t = RR_t!(f3,x,β) #we use lnK here.
        t_inv = t
        t_inv .= 1 ./ t
        df .= 0
        ix = RR_deactivated_phase!(deactivated_phase,β)
        for j in 1:(np-1)
            Kj = viewn(x,nc,j)
            for i in 1:nc
                df[j] += -expm1(Kj[i])*z[i]*t_inv[i]
            end
        end
        RR_deactivate_phase_gh!(ix,df,nothing)
        return df
    end

    function fg(df,β)
        _,_,_,F_cache,f1,f2,f3,dem_cache,Hϕ = ss_cache
        result, idx_vapour = _result
        nc,np = length(z),numphases(result)
        ix = RR_deactivated_phase!(deactivated_phase,β)
        t = RR_t!(f3,x,β) #we use lnK here.
        fx = zero(eltype(β))
        for i in 1:nc
            fx -= z[i]*log(abs(t[i]))
        end
        df .= 0
        t_inv = t
        t_inv .= 1 ./ t
        for j in 1:(np-1)
            Kj = viewn(x,nc,j)
            for i in 1:nc
                df[j] += -expm1(Kj[i])*z[i]*t_inv[i]
            end
        end
        RR_deactivate_phase_gh!(ix,df,nothing)
        return fx,df
    end

    function fgh(df,d2f,β)
        _,_,_,F_cache,f1,f2,f3,dem_cache,Hϕ = ss_cache
        result, idx_vapour = _result
        nc,np = length(z),numphases(result)
        ix = RR_deactivated_phase!(deactivated_phase,β)
        t = RR_t!(f3,x,β) #we use lnK here.
        fx = zero(eltype(β))
        for i in 1:nc
            fx -= z[i]*log(abs(t[i]))
        end
        t_inv = t
        t_inv .= 1 ./ t
        d2f .= 0
        df .= 0
        for j in 1:(np-1)
            Kj = viewn(x,nc,j)
            for i in 1:nc
                ΔKij = expm1(Kj[i])
                ti_inv = t_inv[i]
                zi = z[i]
                df[j] += -ΔKij*zi*ti_inv
                for k in 1:(np-1)
                    Kk = viewn(x,nc,k)
                    ΔKik = expm1(Kk[i])
                    d2f[j,k] += ΔKik*ΔKij*zi*ti_inv*ti_inv
                end
            end
        end
        #recommended by michelsen, to avoid making the hessian singular
        for i in 1:length(β)
            d2f[i,i] += 1e-10
        end
        RR_deactivate_phase_gh!(ix,df,d2f)
        return fx,df,d2f
    end

    function h(d2f,β)
        _,_,_,F_cache,f1,f2,f3,dem_cache,Hϕ = ss_cache
        result, idx_vapour = _result
        nc,np = length(z),numphases(result)
        ix = RR_deactivated_phase!(deactivated_phase,β)
        t = RR_t!(f3,x,β)
        d2f .= 0
        t_inv = t
        t_inv .= 1 ./ t
        for j in 1:(np-1)
            Kj = viewn(x,nc,j)
            for i in 1:nc
                ΔKij = expm1(Kj[i])
                ti_inv = t_inv[i]
                zi = z[i]
                for k in 1:(np-1)
                    Kk = viewn(x,nc,k)
                    ΔKik = expm1(Kk[i])
                    d2f[j,k] += ΔKik*ΔKij*zi*ti_inv*ti_inv
                end
            end
        end
        #recommended by michelsen, to avoid making the hessian singular
        for i in 1:length(β)
            d2f[i,i] += 1e-10
        end
        RR_deactivate_phase_gh!(ix,nothing,d2f)
        return d2f
    end

    return ScalarObjective(f=f,
    g=g,
    fg=fg,
    fgh=fgh,
    h=h)
end

function multiphase_lnK_K!(F,x,np,nc)
    for l in 1:(np - 1)
        lnK = viewn(x,nc,l)
        K = viewn(F,nc,l)
        K .= exp.(lnK)
    end
    return F
end

function multiphase_K_lnK!(F,x,np,nc)
    for l in 1:(np - 1)
        K = viewn(x,nc,l)
        lnK = viewn(F,nc,l)
        lnK .= log.(K)
    end
    return F
end

function multiphase_RR_lnK!(F, x, model, p, T, z, _result, ss_cache)
    _,_,_,F_cache,f1,f2,f3,dem_cache,Hϕ = ss_cache
    result, idx_vapour = _result
    comps, βi, volumes = result.compositions,result.fractions,result.volumes
    nc = length(z)
    np = length(βi)
    imin,βmin = findmin(βi)
    βmin < 0 && return F
    multiphase_lnK_K!(F_cache,x,np,nc) #F_cache now contains Kij
    β = viewn(F,nc,np)

    ##okuno et al. (2010) : calculate t
    t = RR_t!(f3,x,β)
    #calculate new compositions for phase np
    xnp = comps[np]
    xnp .= z ./ t
    xnp ./= sum(xnp)
    phasenp = idx_vapour[] == np ? :vapour : :liquid
    lnϕnp_temp,vnp = modified_lnϕ(model, p, T, xnp, Hϕ, vol0 = last(volumes),phase = phasenp)
    #caches: f3 stores t, f2 stores lnϕ1
    lnϕnp = f2
    lnϕnp .= lnϕnp_temp
    volumes[end] = vnp
    #calculate new compositions for the rest of phases
    for i in 1:(np - 1)
        Ki = viewn(F_cache,nc,i)
        xi = comps[i]
        RR_t_x!(xi,Ki,t,z) #xi .= Ki .* z ./ t, handling the case of K = Inf
        xi ./= sum(xi)
    end

    #update lnK
    for i in 1:np - 1
        lnKi = viewn(F,nc,i)
        xi = comps[i]
        phasei = idx_vapour[] == i ? :vapour : :liquid
        lnϕi,vi = modified_lnϕ(model, p, T, xi, Hϕ, vol0 = volumes[i],phase = phasei)
        volumes[i] = vi
        lnKi .=  lnϕnp .- lnϕi
    end
    return F
end
#multialgorithm to add a new phase.

#if there is only one phase, tries diffusive stability - K values.
#otherwise, tries to split the phases according to tpd and gibbs optimization.
function _add_phases!(model,p,T,z,_result,cache,options)
    result, idx_vapour = _result
    comps, β, volumes = result.compositions,result.fractions,result.volumes
    phase_cache,ss_cache = cache
    tpd_cache,found_tpd,found_tpd_lnphi,found_tpd_volumes = phase_cache
    _,_,_,_,_,Hϕ = tpd_cache
    result_cache = ss_cache[1]
    comp_cache = result_cache[1]
    np = length(comps)
    nc = length(z)
    verbose = options.verbose

    δn_add = false

    np == nc && return false

    max_phases = min(options.max_phases,nc)
    #np >= max_phases && return 0 #we cannot add new phases here
    iter = np
    phases_comps = fill(:unknown,length(comps))
    if idx_vapour[] != 0
        fill!(phases_comps,:liquid)
        phases_comps[idx_vapour[]] = :vapour
    end
    phases_tpd = fill(:unknown,length(found_tpd))
    #step 1: find new tpds.
    for i in 1:iter
        w = comps[i]
        vw = volumes[i]

        if idx_vapour[] == 0
            phase_wi = identify_phase(model,p,T,w,vol = vw)
            if is_vapour(phase_wi)
                idx_vapour[] = i
            end
        else
            if idx_vapour[] == i
                phase_wi = :vapour
            else
                phase_wi = :liquid
            end
        end

        # If we only have a vapour phase and no liquid yet, we must search VLE, not LLE.
        # Restrict to LLE only once we already have both vapour and at least one liquid phase.
        is_lle = (idx_vapour[] != 0) && (length(comps) > 1)
        tpd_i = tpd(model,p,T,w,tpd_cache,reduced = true, strategy = :pure,break_first = true,lle = is_lle)
        # println("[MPFLASH] _add_phases!: i=", i, ", idx_vapour=", idx_vapour[], ", is_lle=", is_lle, ", found=", length(tpd_i[1]))
        length(tpd_i[1]) == 0 && continue
        already_found = false
        tpd_comps,_ttpd,phases_z,phases_w = tpd_i
        if is_vapour(phase_wi) && idx_vapour[] == 0
            idx_vapour[] = i
            fill!(phases_comps,:liquid)
            phases_comps[idx_vapour[]] = :vapour
        end
        y = tpd_comps[1]
        for k in 1:length(found_tpd)
            if z_norm(found_tpd[k],y) < 1e-5
                already_found = true #this tpd result was already found.
            end
        end

        already_found && continue
        phases_comps[i] = phases_z[1]
        push!(phases_tpd,phases_w[1])
        push!(found_tpd,y)
        lnϕy,vy = modified_lnϕ(model,p,T,y,Hϕ,phase = phases_w[1])
        push!(found_tpd_volumes,vy)
        push!(found_tpd_lnphi,copy(lnϕy))
    end
    #step 2: calculate a matrix of tpd for each value.
    # if no TPD candidates were found, nothing to add
    if length(found_tpd) == 0
        verbose && @info("[MPFLASH] _add_phases!: no TPD candidates found")
        return false
    end

    tpds = zeros(length(found_tpd),length(comps))
    #cache di of components
    for i in 1:length(comps)
        di = comp_cache[i]
        __di,_ = modified_lnϕ(model, p, T, comps[i], Hϕ, vol0 = volumes[i])
        di .= __di .+ log.(comps[i])
    end

    for i in 1:length(found_tpd)
        lnϕw = found_tpd_lnphi[i]
        w = found_tpd[i]
        for j in 1:length(comps)
            dj = comp_cache[j]
            tpds[i,j] = 1 + @sum(w[k]*(log(w[k]) + lnϕw[k] - dj[k] - 1))
        end
    end
    #step 3: find the tpd with the lowest tpd value
    # treat extremely small negative tpd as zero to avoid degenerate splits
    tpd_tol = 1e-12
    # println("[MPFLASH] _add_phases!: min(tpds)=", minimum(tpds), ", tol=", tpd_tol)
    (minimum(tpds) > -tpd_tol) && return false
    #check our current phases and the trial ones
    for ii in 1:length(found_tpd)
        for jj in 1:length(comps)
            tpds[ii,jj] > 0 && continue
            idx_vapour[] == jj && length(comps) == 2 && continue #only check liquid when we have a liquid and a vapour
            w,vw = comps[jj],volumes[jj]

            #identify phase if not done
            if is_unknown(phases_comps[jj]) && idx_vapour[] == 0
                phases_comps[jj] = identify_phase(model,p,T,w,vol = vw)
            end

            phase_w = phases_comps[jj]
            if is_vapour(phase_w) && idx_vapour[] == 0 #we identified the vapour phase
                idx_vapour[] = jj
            end

            y,vy = found_tpd[ii],found_tpd_volumes[ii]
            if is_unknown(phases_tpd[ii])
                phases_tpd[ii] = identify_phase(model,p,T,y,vol = vy)
            end
            phase_y = phases_tpd[ii]
            #vy = volume(model,p,T,y,phase = phase_y)
            #phase not stable: generate a new one from tpd result
            result_split_tpd = split_phase_tpd(model,p,T,w,y,phase_w,phase_y,vw,vy)
            β1,β2 = result_split_tpd.fractions
            x1,x2 = result_split_tpd.compositions
            v1,v2 = result_split_tpd.volumes
            dgi = result_split_tpd.data.g
            # println("[MPFLASH] split result: β1=",β1,", β2=",β2,", dgi=",dgi)
            # println("[MPFLASH] split x1 sum=",sum(x1),", any_nan=",any(!isfinite,x1),", v1=",v1)
            # println("[MPFLASH] split x2 sum=",sum(x2),", any_nan=",any(!isfinite,x2),", v2=",v2)
            #check that the new generated phase is not equal to one existing composition
            knew = 0
            for i in 1:length(comps)
                if z_norm(comps[i],x1) < 1e-5 && abs(1/v1 - 1/volumes[i]) <= 1e-4
                    knew = i #the split resulted in a new phase equal to one already existing
                end
            end

            if (!isnan(dgi) && (dgi < 0) && !iszero(β2)) || isone(np) #new phase found
                β0 = β[jj]
                β[jj] = β0*β2
                comps[jj] = x2
                volumes[jj] = v2
                if iszero(knew)
                    #new incipient phase
                    push!(comps,x1)
                    push!(volumes,v1)
                    push!(β,β0*β1)
                    if is_vapour(identify_phase(model,p,T,x1,vol = v1))
                        idx_vapour[] = length(comps)
                    end
                    return true
                else
                    #the new phase is similar to an existing one, we redistribute the fractions and volumes, without adding new phases.
                    wi,wj = x1,comps[jj]
                    vi,vj = v1,volumes[jj]
                    βi,βj = β0*β1,β[jj]
                    volumes[jj] = (βi*vi + βj*vj)/(βi + βj)
                    β[jj] = (βi + βj)
                    wj .= (βi .* wi .+ βj .* wj) ./ (βi .+ βj)
                    return true
                end
            end
        end
    end
    return false
end

#find the phase with the minimum βi. if mixing that phase with any other phase generates
#a more stable phase, remove it
function _remove_phases!(model,p,T,z,_result,cache,options)
    result, idx_vapour = _result
    comps, β, volumes = result.compositions,result.fractions,result.volumes
    δn_remove = false
    phase_cache,ss_cache = cache
    tpd_cache,found_tpd,found_tpd_lnphi,found_tpd_volumes = phase_cache
    #strategy 0: remove all "equal" phases.
    #if phases are equal (equal volume and comps), fuse them

    np = numphases(result)
    for i in 1:(np*np)
        i,j = findfirst_duplicate_phases(result,false)
        if i == j == 0
            break
        end
        merge_phase!(result,i,j)
        if idx_vapour[] == j
            idx_vapour[] = i
        end
    end

    if numphases(result) == 1
        return true
    end

    #strategy A: remove all phases with βi < βmin = 4eps(eltype(βi))
    β_remove = findall(<(4eps(eltype(β))),β)
    adjust_idx_vapour!(idx_vapour,β_remove)

    if length(β_remove) > 0
        #remove all phases with negative values
        delete_phase!(result,β_remove)

        #add a new phase and reconstitute it
        push!(comps,similar(comps[1]))
        push!(β,0.0)

        reconstitute_x!(comps,z,β,length(β))
        push!(volumes,volume(model,p,T,comps[end]))
        if idx_vapour[] == 0
            if is_vapour(identify_phase(model,p,T,comps[end],vol = volumes[end]))
                idx_vapour[] = length(comps)
            end
        end
        return true
    end

    #strategy B: remove one phase that does not help in equilibria
    βmin,imin = findmin(β)
    wmin,vmin = comps[imin],volumes[imin]
    phase_min = __mpflash_phase(idx_vapour[],imin)
    gmin,_ = modified_gibbs(model,p,T,wmin,phase_min,vmin)
    wmix = similar(comps[1])
    for i in 1:np
        i == imin && continue
        wi,βi,vi = comps[i],β[i],volumes[i]
        wmix .= βi .* wi .+ βmin .* wmin
        βmix = βi + βmin
        wmix .= wmix ./ sum(wmix)
        phasei = __mpflash_phase(idx_vapour[],i)
        #gi = eos(model,vi,T,wi) + vi*p
        gi,_ = modified_gibbs(model,p,T,wi,phasei,vi)
        #gmix = eos(model,vmix,T,wmix) + vmix*p
        gmix,vmix = modified_gibbs(model,p,T,wmix,:unknown)
        Δg = βmix*gmix - βi*gi - βmin*gmin
        if Δg < 0 #the mixed phase has a lower Gibbs energy than the sum of its parts. remove minimum fraction phase.
            δn_remove = true
            comps[i] = wmix
            β[i] = βmix
            volumes[i] = vmix
            delete_phase!(result,imin)
            adjust_idx_vapour!(idx_vapour,imin)
            break
        end
    end
    return false
end

#given a list of β_remove, remove the corresponding phase
function adjust_idx_vapour!(idx_vapour,β_remove)
    if idx_vapour[] in β_remove
        idx_vapour[] = 0
    end

    idx_vapour_new = idx_vapour[]
    idx_vapour_current = idx_vapour[]
    #when removing phases, adjust the index of the vapour phase
    if idx_vapour_current != 0
        for β_idx in β_remove
            if idx_vapour_current >= β_idx
                idx_vapour_new = idx_vapour_new - 1
            end
        end
        idx_vapour[] = idx_vapour_new
    end
    return idx_vapour[]
end

function set_idx_vapour!(idx_vapour,model,result)
    idx_vapour[] != 0 && return idx_vapour[]

    np = numphases(result)
    p,T = result.data.p,result.data.T
    v = result.volumes
    β = result.fractions
    x = result.compositions

    for i in 1:np
        if is_vapour(identify_phase(model,p,T,x[i],vol = v[i]))
            idx_vapour[] = i
            break
        end
    end
    return idx_vapour[]
end



function split_phase_tpd(model,p,T,z,w,phase_z = :unknown,phase_w = :unknown,vz = volume(model,p,T,z,phase = phase_z),vw = volume(model,p,T,w,phase = phase_w))
    #w is a phase found via tpd.
    #=if w has a negative tpd, then it exists phases x1,x2 such as
    β*g(x1) + (1 - β)*g(x2) < g(w)
    x1 = w
    x2 = (z -  β*w)/(1 - β)


    βmin = 0
    z -  β*w = 0
    z/w = β
    βmax = minimum(zi/wi for i in comps)
    =#
    kk = (z = z,w = w)
    #@show kk
    β1 = zero(eltype(w))
    β2 = one(eltype(w))
    # Robust upper bound: consider only indices with w[i] > 0 to avoid 0/0 and Inf pollution
    for (zi, wi) in zip(z,w)
        β2 = ifelse(wi > 0, min(β2, zi/wi), β2)
    end
    x2 = similar(w)
    g1,_ = modified_gibbs(model,vw,T,w,phase_w,vw)
    gz,_ = modified_gibbs(model,vw,T,z,phase_z,vz)

    if isone(β2)
        x2 .= z
    else
        @. x2 = (z -  β2 * w) / (1 - β2)
    end
    #@show β2
    x3 = x2
    phase = phase_w == phase_z ? phase_z : (is_vapour(phase_w) ? :liquid : :unknown)
    cache = Ref(g1)

    function ff(_β)
        x3 .= (z .-  _β .* w) ./ (1 .- _β)
        gw,_vw = modified_gibbs(model,p,T,x3,phase)
        dgi = _β*g1 + (1-_β)*gw - gz
        cache[] = _vw
        return dgi
    end

    sol = Solvers.optimize(ff,(zero(β2),β2),Solvers.BoundOptim1Var(),NLSolvers.OptimizationOptions(x_abstol = 1e-5))
    βi = Solvers.x_sol(sol)
    dgi_sol = Solvers.x_minimum(sol)
    v3_sol = cache[]
    z2 = z
    x3 .= (z .-  βi .* w) ./ (1 .- βi)
    fracs = SVector((1 - βi),βi)
    volumes = SVector(v3_sol,vw)
    comps = SVector(x3,w)
    return FlashResult(comps,fracs,volumes,FlashData(p,T,dgi_sol))
end


function tp_flash_multi_neq!(model,p,T,z,_result,ss_cache,options)
    result_cache,x0,x,xx,f1,f2,f3,dem_cache,Hϕ = ss_cache
    result, idx_vapour = _result
    comps, βi, volumes = result.compositions,result.fractions,result.volumes
    np = numphases(result)
    nc = length(z)
    resize!(x,(np-1)*nc + (np - 1) + np)
    vx = @view x[(end - np + 1):end]
    for i in 1:np
        vx[i] = log(volumes[i])
    end
    if any(!isfinite,βi) || minimum(βi) < 0
        gmix = zero(eltype(βi))/0
        false,gmix
    end
    opt_options = OptimizationOptions(maxiter = 30)
    f = multi_g_obj(model,p,T,z,_result,ss_cache)
    if options.second_order
        sol = Solvers.optimize(f,x,LineSearch(Newton2(x),Backtracking()),opt_options)
        gmix = Solvers.x_minimum(sol)
        F = Solvers.x_sol(sol)
    else
        sol = Solvers.optimize(f,x,LineSearch(LBFGS()))
        gmix = Solvers.x_minimum(sol)
        F = Solvers.x_sol(sol)
    end

    idx_β_begin = (np-1)*nc + 1
    idx_β_end = idx_β_begin + np - 2
    β = view(F,idx_β_begin:idx_β_end)
    v = @view F[(end - np + 1):end]
    xi = f1
    xnp = f2
    t = RR_t!(f1,F,β)
    βnp = 1 - sum(β)
    xnp .= z ./ t
    xx_np = comps[np]
    xx_np .= xnp
    volumes[np] = exp(v[np])
    for i in 1:np-1
        xxi = comps[i]
        Ki = viewn(F,nc,i)
        xxi .= xnp .* exp.(Ki)
        βi[i] = β[i]
        volumes[i] = exp(v[i])
    end
    βi[np] = βnp
    converged = neq_converged(model,p,T,z,result)
    return converged,gmix
end

function multi_g_obj(model,p,T,z,_result,ss_cache)
    result_cache,x0,x,xx,f1,f2,f3,dem_cache,Hϕ = ss_cache
    result, idx_vapour = _result
    nc = length(z)
    np = numphases(result)
    function f(𝕏)
        xnp = similar(𝕏,nc)
        xi = similar(𝕏,nc)
        idx_β_begin = (np-1)*nc + 1
        idx_β_end = idx_β_begin + np - 2
        β = view(𝕏,idx_β_begin:idx_β_end)
        vols = view(𝕏,(idx_β_end+1):length(𝕏))
        t = RR_t!(xi,𝕏,β)
        βnp = 1 - sum(β)
        xnp .= z ./ t
        vnp = vols[np]
        phase_np = __mpflash_phase(idx_vapour[],np)

        #g = βnp*(eos(model,vnp,T,xnp) + p*vnp)
        if has_a_res(model)
            g = βnp*(eos(model,vnp,T,xnp) + p*vnp)
        else
            g = βnp*modified_gibbs(model,p,T,xnp,phase_np,vnp)[1]
        end

        for i in 1:np-1
            Ki = viewn(𝕏,nc,i)
            xi .= xnp .* exp.(Ki)
            vi = exp(vols[i])
            phase_i = __mpflash_phase(idx_vapour[],i)
            #g += β[i]*(eos(model,vi,T,xi) + p*vi)
            if has_a_res(model)
                g += β[i]*(eos(model,vi,T,xi) + p*vi)
            else
                g += β[i]*modified_gibbs(model,p,T,xi,phase_i,vi)[1]
            end
        end
        return g/(Rgas(model)*T)
    end

    return f
end

function initial_beta!(comps,z)
    isempty(comps) && throw(ArgumentError("initial_beta! received no phase compositions (comps is empty)"))
    βi = reduce(hcat,comps) \ z
    βmin,imin = findmin(βi)

    if βmin < 0
        βi[imin] = 0
        reconstitute_x!(comps,z,βi,imin)
    end
    βi ./= sum(βi)
    return βi
end

function reconstitute_x!(comps,z,bi,i0)
    #we suppose that bi is already set to zero.
    z1 = similar(comps[1])
    z1 .= 0
    for i in 1:length(bi)
        if i != i0
            z1 .+= bi[i] * comps[i]
        end
    end
    z1 ./= sum(z1)
    β = one(eltype(z1))
    z2 = comps[i0]
    for i in 1:length(z)
        β = min(β,z[i]/z1[i])
    end
    β = 0.5*β
    z2 .= (z .- β .*z1)/(1 .- β)
    bi .*= (1 - β)
    bi[i0] = β
    #we split the phases between z and z1.
end