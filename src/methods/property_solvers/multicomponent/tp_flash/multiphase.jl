struct MultiPhaseTPFlash{T} <: TPFlashMethod
    K0::Union{Vector{T},Nothing}
    K_tol::Float64
    ss_iters::Int
    nacc::Int
    second_order::Bool
    max_phases::Int
    phase_iters::Int
end

"""
    MultiPhaseTPFlash(;kwargs...)

Method to solve non-reactive multiphase (`np` phases), multicomponent (`nc` components) flash problem.

The flash algorithm uses successive stability tests to find new phases [1], and then tries to solve the system via rachford-rice and succesive substitution for `nc * np * ss_iters` iterations.

If the Rachford-Rice SS fails to converge, it proceeds to solve the system via gibbs minimization in VT-space using lnK-β-ρ as variables [3].

The algorithm finishes when SS or the gibbs minimization converges and all resulting phases are stable.

If the result of the phase equilibria is not stable, then it proceeds to add/remove phases again, for a maximum of `phase_iters` iterations.

### Keyword Arguments:

- `K0` (optional), initial guess for the constants K
- `x0` (optional), initial guess for the composition of phase x
- `y0` (optional), initial guess for the composition of phase y
- `K_tol = sqrt(eps(Float64))`, tolerance to stop the calculation (`norm(lnK,1) < K_tol`)
- `ss_iters = 4`, number of Successive Substitution iterations to perform
- `nacc = 3`, accelerate successive substitution method every nacc steps. Should be a integer bigger than 3. Set to 0 for no acceleration.
- `second_order = true`, whether to solve the gibbs energy minimization using the analytical hessian or not. If set to `false`, the gibbs minimization will be done using L-BFGS.
- `max_phases = typemax(Int)`, the algorithm stops if there are more than `min(max_phases,nc)` phases
- `phase_iters = 20`, the maximum number of solve-add/remove-phase iterations

## References
1. Thermopack - Thermodynamic equilibrium algorithms reimplemented in a new framework. (2020, September 08). https://github.com/thermotools/thermopack. Retrieved May 4, 2024, from https://github.com/thermotools/thermopack/blob/main/docs/memo/flash/flash.pdf
2.  Okuno, R., Johns, R. T. T., & Sepehrnoori, K. (2010). A new algorithm for Rachford-Rice for multiphase compositional simulation. SPE Journal, 15(02), 313–325. [doi:10.2118/117752-pa](https://doi.org/10.2118/117752-pa)
3.  Adhithya, T. B., & Venkatarathnam, G. (2021). New pressure and density based methods for isothermal-isobaric flash calculations. Fluid Phase Equilibria, 537(112980), 112980. [doi:10.1016/j.fluid.2021.112980](https://doi.org/10.1016/j.fluid.2021.112980)
"""
function MultiPhaseTPFlash(;
    K0 = nothing,
    x0 = nothing,
    y0 = nothing,
    K_tol = sqrt(eps(Float64)),
    ss_iters = 2,
    nacc = 5,
    second_order = true,
    max_phases = typemax(Int),
    phase_iters = 20) #TODO: find a better value for this)
    if K0 == x0 == y0 == nothing #nothing specified
    #is_lle(equilibrium)
        T = Nothing
    else
        if !isnothing(K0) & isnothing(x0) & isnothing(y0) #K0 specified
            T = eltype(K0)
        elseif isnothing(K0) & !isnothing(x0) & !isnothing(y0)  #x0, y0 specified
            K0 = x0 ./ y0
            T = eltype(K0)
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

    return MultiPhaseTPFlash{T}(K0,K_tol,ss_iters,nacc,second_order,max_phases,phase_iters)
end

export MultiPhaseTPFlash

function index_reduction(m::MultiPhaseTPFlash,idx::AbstractVector)
    K0,K_tol,ss_iters,nacc,second_order,max_phases,phase_iters = m.K0,m.K_tol,m.ss_iters,m.nacc,m.second_order,m.max_phases,m.phase_iters
    K0 !== nothing && (K0 = K0[idx])
    T = eltype(K0)
    return MultiPhaseTPFlash{T}(K0,K_tol,ss_iters,nacc,second_order,max_phases,phase_iters)
end

function tpd_cache end

function tp_flash_multi_cache(model,p,T,z)
    pure = split_model(model)
    _tpd_cache = tpd_cache(model,p,T,z)

    F_cache = zeros(Base.promote_eltype(model,p,T,z),length(z))
    F_cache2 = similar(F_cache)
    F_cache3 = similar(F_cache)
    f1,f2,f3 = similar(F_cache),similar(F_cache),similar(F_cache)
    F3,F4,F5,ΔF1,ΔF2,xdem,Fdem = similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache)
    dem_cache = (F3,F4,F5,ΔF1,ΔF2,xdem,Fdem)
    comps_cache = fill(similar(f1),1)
    v_cache = similar(f1)
    bi_cache = similar(f1)
    result_cache = comps_cache,bi_cache,v_cache
    lnφv,lnφl = similar(f1),similar(f1)
    for i in 1:length(z)
        vl = volume(pure[i],p,T,phase = :l)
        vv = volume(pure[i],p,T,phase = :v)
        lnφv[i] = VT_lnϕ_pure(pure[i],vv,T,p)
        lnφl[i] = VT_lnϕ_pure(pure[i],vl,T,p)
    end
    phase_cache = (pure,_tpd_cache,lnφv,lnφl)
    ss_cache = result_cache,F_cache,F_cache2,F_cache3,f1,f2,f3,dem_cache
    return phase_cache,ss_cache,nothing
end

function resize_cache!(cache,np)
    phase_cache,ss_cache,newton_cache = cache
    nc = length(phase_cache[1])
    result_cache,F_cache,F_cache2,F_cache3,f1,f2,f3,dem_cache = ss_cache
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
    return tp_flash_multi(model,p,T,z,method)
end

function tp_flash_multi(model,p,T,z,options = MultiPhaseTPFlash())
    n_phases = 1
    max_iter = options.phase_iters
    V = p
    z̄ = zeros(@f(Base.promote_eltype),length(z))
    z̄ .= z
    z̄ ./= sum(z)
    comps = [z̄]
    βi = [one(eltype(z))]
    vz = volume(model,p,T,z)
    volumes = [vz]

    idx_vapour = Ref(0)
    _result = (comps, βi, volumes, idx_vapour)
    result = (comps, βi, volumes)
    δn_remove = 0

    cache = tp_flash_multi_cache(model,p,T,z)
    phase_cache,ss_cache,newton_cache = cache
    done = false
    converged = false
    iter = 0
    ss_converged = false
    neq_converged = false
    #step 1: tpd: we check if there is any unstable phase.
    g0 = (eos(model,volumes[1],T,z) + p*vz)/(Rgas(model)*T)
    δn_add = _add_phases!(model,p,T,z,_result,cache,options)
    if δn_add == 0
        return comps, βi, volumes,g0
    end
    gmix = NaN*g0
    #step 2: main loop,iterate flashes until all phases are stable
    while !done
        iter += 1
        δn = δn_add - δn_remove
        if δn_add != 0 || δn_remove != 0#step 2.1: δn != 0. we add or remove phases and find new candidate ones.
            δn_add,δn_remove = 0,0
            n_phases = length(comps)
            #sucessive substitution iteration
            cache = resize_cache!(cache,n_phases)
            _result,ss_converged = tp_flash_multi_ss!(model,p,T,z,_result,ss_cache,options)
            if !ss_converged
                _result,neq_converged,gmix = tp_flash_multi_neq!(model,p,T,z,_result,ss_cache,options)
            end

            δn_add = _add_phases!(model,p,T,z,_result,cache,options)
            if δn_add != 1
                δn_remove = _remove_phases!(model,p,T,z,_result,cache,options)
            end

            converged = neq_converged || ss_converged
            no_new_phases = δn_add == δn_remove == 0
            converged = converged && no_new_phases

            converged && ss_converged && (gmix = _multiphase_gibbs(model,p,T,result)/(Rgas(model)*T))
        end
        done = iter > max_iter
        done = done || converged
    end

    return comps, βi, volumes, gmix
end

function neq_converged(model,p,T,z,result)
    #TODO: add convergence criteria
    return true
end

function tp_flash_multi_ss!(model,p,T,z,_result,ss_cache,options)
    result_cache,x0,x,xx,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes, idx_vapour = _result
    result = comps, βi, volumes
    _result_cache = result_cache[1],result_cache[2],result_cache[3],idx_vapour
    nc = length(z)
    np = length(βi)
    #set the initial values to 0
    x0 .= 0
    x .= 0
    xx .= 0
    comps_cache,βi_cache,volumes_cache = result_cache
    #fill F with data
    xnp = comps[np]
    β = viewn(x0,nc,np)
    β = view(β,1:np)
    β .= βi
    for i in 1:(np-1)
        lnKi = viewn(x0,nc,i)
        xi = comps[i]
        lnKi .= log.(xi ./ xnp)
    end

    #options
    ss_iters = options.ss_iters
    K_tol = options.K_tol
    nacc = options.nacc

    max_iters = ss_iters*np*nc
    itacc = 0
    converged = false
    F3,F4,F5,ΔF1,ΔF2,xdem,Fdem = dem_cache
    for i in 1:max_iters
        itacc += 1
        fixpoint_multiphase!(x, x0, model, p, T, z, _result, ss_cache)
        #store values in cache
        βi_cache .= βi
        volumes_cache .= volumes
        for i in 1:np
            xi_cache,xi = comps_cache[i],comps[i]
            xi_cache .= xi
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
            fixpoint_multiphase!(Fdem, xdem, model, p, T, z, _result_cache, ss_cache)

            gibbs = _multiphase_gibbs(model,p,T,result)
            gibbs_dem = _multiphase_gibbs(model,p,T,result_cache)
            if gibbs_dem < gibbs
                βi .= βi_cache
                volumes .= volumes_cache
                for i in 1:np
                    xi,xi_cache = comps[i],comps_cache[i]
                    xi .= xi_cache
                end
                x .= Fdem
            end
        end
        converged = dnorm(x0,x,1) < K_tol
        converged && break
        x0 .= x
    end

    return _result,converged
end

function fixpoint_multiphase!(F, x, model, p, T, z, result, ss_cache)
    #given constant K, calculate β
    multiphase_RR_β!(F, x, z, result, ss_cache)
    #given constant β, calculate K
    multiphase_RR_lnK!(F, x, model, p, T, z, result, ss_cache)
    return F
end

function RR_t!(t,x,β,np,nc)
    t .= 1
    for l in 1:(np - 1)
        Kl = viewn(x,nc,l)
        βl = β[l]
        fi = zero(βl)
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

#constant K, calculate β
function multiphase_RR_β!(F, x, z, _result, ss_cache)
    _,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes, idx_vapour = _result
    result = comps, βi, volumes
    nc = length(z)
    np = length(βi)
    βF = viewn(F,nc,np)

    #transforms lnK to K
    multiphase_lnK_K!(F_cache,x,np,nc)
    if np == 2
        K = viewn(F_cache,nc,1)
        βsol = rachfordrice(K,z)
        βF[1] = βsol #two-phase solution does not require solving a neq problem
        βF[2] = 1 - βsol
        βi[1] = βsol
        βi[2] = 1 - βsol
        return F
    end

    f = RR_obj(x, z, _result, ss_cache)
    β0 = viewn(x,nc,np)
    β0 = copy(view(β0,1:np-1))
    #=
    b = f2
    b .= 1 .- z
    for j in 1:(np-1)
        Kj = viewn(F_cache,nc,j)
        for i in 1:nc
            b[i] = min(b[i],1 - Kj[i]*z[i])
        end
    end
    a = view(f1,1:np-1)
    #a = similar(b,np-1)

    function ls_restricted(φ::P,λ) where P
        d = φ.d
        x = φ.z
        λmax = λ
        for i in 1:nc
            ai = RR_β_ai!(a,F_cache,i,nc,np)
            ai_d = dot(ai,d)
            if ai_d > 0
                λx = one(λ)*(b[i] - dot(ai,x))/ai_d
                λmax = min(λx,λmax)
            end
        end
        return min(λmax,one(λmax))
    end =#
    ls_restricted(x1,x2) = x2
    #ls = RestrictedLineSearch(ls_restricted,Static(1.0))
    βresult = Solvers.optimize(f,β0,LineSearch(Newton(linsolve = static_linsolve)))
    βsol = Solvers.x_sol(βresult)
    ∑β = sum(βsol)
    for i in 1:(np-1)
        βF[i] = βsol[i]
        βi[i] = βsol[i]
    end
    #@show βsol
    #update last phase fraction
    βF[np] = 1 - ∑β
    βi[np] = 1 - ∑β
    return F
end

function RR_obj(x, z, _result, ss_cache)
    #okuno et al. (2010): objetive function, gradient, hessian
    function f(β)
        _,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
        comps, βi, volumes, idx_vapour = _result
        nc,np = length(z),length(comps)
        t = RR_t!(f3,x,β,np,nc) #we use lnK here.
        fx = zero(eltype(β))
        for i in 1:nc
            fx -= z[i]*log(abs(t[i]))
        end
        return fx
    end

    function g(df,β)
        _,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
        comps, βi, volumes, idx_vapour = _result
        nc,np = length(z),length(comps)
        t = RR_t!(f3,x,β,np,nc) #we use lnK here.
        t_inv = t
        t_inv .= 1 ./ t
        df .= 0
        for j in 1:(np-1)
            Kj = viewn(x,nc,j)
            for i in 1:nc
                df[j] += -expm1(Kj[i])*z[i]*t_inv[i]
            end
        end
        return df
    end

    function fg(df,β)
        _,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
        comps, βi, volumes, idx_vapour = _result
        nc,np = length(z),length(comps)
        t = RR_t!(f3,x,β,np,nc) #we use lnK here.
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
        return fx,df
    end

    function fgh(df,d2f,β)
        _,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
        comps, βi, volumes, idx_vapour = _result
        nc,np = length(z),length(comps)
        t = RR_t!(f3,x,β,np,nc) #we use lnK here.
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
        return fx,df,d2f
    end

    function h(d2f,β)
        _,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
        comps, βi, volumes, idx_vapour = _result
        nc,np = length(z),length(comps)
        t = RR_t!(f3,x,β,np,nc)
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
    _,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes, idx_vapour = _result
    result = comps, βi, volumes
    nc = length(z)
    np = length(βi)
    multiphase_lnK_K!(F_cache,x,np,nc) #F_cache now contains Kij
    β = viewn(F,nc,np)

    ##okuno et al. (2010) : calculate t
    t = RR_t!(f1,x,β,np,nc)

    #calculate new compositions for phase np
    xnp = comps[np]
    xnp .= z ./ t
    xnp ./= sum(xnp)
    lnϕnp,vnp = lnϕ!(f2, model, p, T, xnp, vol0 = last(volumes))
    isnan(vnp) && (lnϕnp,vnp = lnϕ!(f2, model, p, T, xnp)) #restart
    volumes[end] = vnp

    #calculate new compositions for the rest of phases
    for i in 1:(np - 1)
        Ki = viewn(F_cache,nc,i)
        xi = comps[i]
        xi .= Ki .* z ./ t
        xi ./= sum(xi)
    end

    #update lnK
    for i in 1:np - 1
        lnKi = viewn(F,nc,i)
        xi = comps[i]
        lnϕi,vi = lnϕ!(f3, model, p, T, xi, vol0 = volumes[i])
        if isnan(vi) #restart
            lnϕi,vi = lnϕ!(f1, model, p, T, xi)
        end
        volumes[i] = vi
        lnKi .=  lnϕnp .- lnϕi
    end
    return F
end
#multialgorithm to add a new phase.

#if there is only one phase, tries diffusive stability - K values.
#otherwise, tries to split the phases according to tpd and gibbs optimization.
function _add_phases!(model,p,T,z,result,cache,options)
    comps, β, volumes,_idx_vapour = result
    idx_vapour = _idx_vapour[]
    phase_cache,ss_cache,newton_cache = cache
    pures,tpd_cache,lnφv,lnφl = phase_cache
    np = length(comps)
    nc = length(z)
    δn_add = 0
    max_phases = min(options.max_phases,nc)
    max_phases,np
    np >= max_phases && return 0 #we cannot add new phases here
    iter = np

    for i in 1:iter
        #np + δn_add == max_phases && break
        w = comps[i]
        vw = volumes[i]
        is_lle = idx_vapour != i && idx_vapour != 0 #if we have a vapour phase, we only search for liquid-liquid splits
        tpd_i = tpd(model,p,T,w,tpd_cache,reduced = true,break_first = true, strategy = :pure,lle = is_lle)

        if length(tpd_i[1]) > 0
            δn_add += 1
            np += 1
            tpd_comps,_,phases_z,phases_w = tpd_i
            #check our current phases and the trial ones
            y = tpd_comps[1]
            phase_w = phases_z[1]
            phase_y = phases_w[1]
            vy = volume(model,p,T,y,phase = phase_y)
            #phase not stable: generate a new one from tpd result
            β1,x1,v1,β2,x2,v2 = split_phase_tpd(model,p,T,w,y,phase_w,phase_y,vw)
            if is_vapour(phase_w) && idx_vapour == 0 #we identified the vapour phase
                _idx_vapour[] = i
            end
            β0 = β[i]
            β[i] = β0*β2
            comps[i] = x2
            volumes[i] = v2
            push!(comps,x1)
            push!(volumes,v1)
            push!(β,β0*β1)
            stable_new = VT_diffusive_stability(model,v1,T,x1)
            if !stable_new
                #the new generated phase is not stable
                #we want to split that
                iter += 1
            end
        end
        δn_add == 1 && break
    end
    return δn_add
end

#find the phase with the minimum βi. if mixing that phase with any other phase generates
#a more stable phase, remove it
function _remove_phases!(model,p,T,z,result,cache,options)
    comps, β, volumes, idx_vapour = result
    n = length(comps)
    δn_remove = 0
    phase_cache,ss_cache,newton_cache = cache
    pures,tpd_cache,lnφv,lnφl = phase_cache
    βmin,imin = findmin(β)
    wmin,vmin = comps[imin],volumes[imin]
    gmin = VT_gibbs_free_energy(model,vmin,T,wmin)
    wmix = similar(comps[1])
    for i in 1:n
        i == imin && continue
        wi,βi,vi = comps[i],β[i],volumes[i]
        wmix .= βi .* wi .+ βmin .* wmin
        βsum = βi + βmin
        wmix .= wmix ./ βsum
        vmix0 = (βi*vi + βmin * vmin)/βsum
        vmix = volume(model,p,T,wmix,vol0 = vmix0)
        gi = eos(model,vi,T,wi) + vi*p
        gmix = eos(model,vmix,T,wmix) + vmix*p
        Δg = βsum*gmix - βi*gi - βmin*gmin
        if Δg < 0 #the mixed phase has a lower gibbs energy than the sum of its parts. remove minimum fraction phase.
            δn_remove += 1
            comps[i] = wmix
            β[i] = βsum
            volumes[i] = vmix
            deleteat!(comps,imin)
            deleteat!(β,imin)
            deleteat!(volumes,imin)
            break
        end
    end
    return δn_remove
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
    β1 = zero(eltype(w))
    β2 = one(eltype(w))

    for i in 1:length(z)
        β2 = min(β2,z[i]/w[i])
    end
    x2 = similar(w)
    #f1 = eos(model,vz,T,z)
    #β = range(β1,β2,length = 20)
    g1 = eos(model,vw,T,w) + p*vw
    gz = eos(model,vz,T,z) + p*vz
    g3 = Inf*g1
    x2 .= (z .-  β2 .* w) ./ (1 .- β2)
    x3 = x2
    phase = phase_w == phase_z ? phase_z : :unknown
    v2 = volume(model,p,T,x2,threaded = false,phase = phase)
    v3 = zero(v2)
    g2 = eos(model,v2,T,x2) + p*v2
    dg1 = g1 - gz
    dg2 = β2*g1 + (1-β2)*g2 - gz
    ϕ = 0.6180339887498949
    βi = ϕ*β1 + (1-ϕ)*β2
    βi0 = one(βi)*10
    _1 = one(βi)

    for i in 1:20
        x3 .= (z .-  βi .* w) ./ (1 .- βi)
        #@show x2
        v3 = volume(model,p,T,x2,threaded = false,phase = phase)
        g3 = eos(model,v2,T,x3) + p*v3
        dgi = βi*g1 + (1-βi)*g3 - gz

        #quadratic interpolation
        A = @SMatrix [β1*β1 β1 _1; β2*β2 β2 _1; βi*βi βi _1]
        b = SVector(dg1,dg2,dgi)
        c = A\b
        βi_intrp = -0.5*c[2]/c[1]
        if dgi < dg2 < dg1
            dg1,β1 = dgi,βi
        elseif dgi < dg1 < dg2
            dg2,β2,v2 = dgi,βi,v3
        elseif dgi < dg1
            dg1,β1,v1 = dgi,βi,v3
        elseif dgi < dg2
            dg2,β2,v2 = dgi,βi,v3
        else
            break
        end
        βi0 = βi
        βi_bisec = ϕ*β1 + (1-ϕ)*β2
        βi = β1 <= βi_intrp <= β2 ? βi_intrp : βi_bisec
        abs(βi0 - βi) < 1e-5 && break
    end
    #@assert βi*w + (1-βi)*x3 ≈ z
    return (1-βi),x3,v3,βi,w,vw
end

function split_phase_k(model,p,T,z,K = nothing,vz = volume(model,p,T,z),pures = split_model(model))
    #split phase via K values.
    if isnothing(K)
        Ki = suggest_K(model,p,T,z,pures)
    else
        Ki = K
    end
    βi = rachfordrice(Ki,z)
    wx = similar(z)
    wy = similar(z)
    rr_flash_liquid!(wx,Ki,z,βi)
    wy .= wx .* Ki
    wx ./= sum(wx)
    wy ./= sum(wy)
    vx = volume(model,p,T,wx)
    vy = volume(model,p,T,wy)
    #@assert βi*wy + (1-βi)*wx ≈ z
    return 1-βi,wx,vx,βi,wy,vy
end

function tp_flash_multi_neq!(model,p,T,z,result,ss_cache,options)
    result_cache,x0,x,xx,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes, idx_vapour = result
    np = length(comps)
    nc = length(z)
    resize!(x,(np-1)*nc + (np - 1) + np)
    vx = @view x[(end - np + 1):end]
    for i in 1:np
        vx[i] = volumes[i]
    end
    f = multi_g_obj(model,p,T,z,result,ss_cache)
    if options.second_order
        sol = Solvers.optimize(f,x,LineSearch(Newton(linsolve = static_linsolve),Backtracking()))
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
    t = RR_t!(f1,F,β,np,nc)
    βnp = 1 - sum(β)
    xnp .= z ./ t
    xx_np = comps[np]
    xx_np .= xnp
    volumes[np] = v[np]
    for i in 1:np-1
        xxi = comps[i]
        Ki = viewn(F,nc,i)
        xxi .= xnp .* exp.(Ki)
        βi[i] = β[i]
        volumes[i] = v[i]
    end
    converged = neq_converged(model,p,T,z,result)
    return result,converged,gmix
end

function multi_g_obj(model,p,T,z,_result,ss_cache)
    result_cache,x0,x,xx,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes, idx_vapour = _result
    nc = length(z)
    np = length(comps)
    function f(𝕏)
        xnp = similar(𝕏,nc)
        xi = similar(𝕏,nc)
        idx_β_begin = (np-1)*nc + 1
        idx_β_end = idx_β_begin + np - 2
        β = view(𝕏,idx_β_begin:idx_β_end)
        vols = view(𝕏,(idx_β_end+1):length(𝕏))
        t = RR_t!(xi,𝕏,β,np,nc)
        βnp = 1 - sum(β)
        xnp .= z ./ t
        vnp = vols[np]
        g = βnp*(eos(model,vnp,T,xnp) + p*vnp)
        for i in 1:np-1
            Ki = viewn(𝕏,nc,i)
            xi .= xnp .* exp.(Ki)
            vi = vols[i]
            g += β[i]*(eos(model,vi,T,xi) + p*vi)
        end
        return g/(Rgas(model)*T)
    end

    return f
end


function split_phase_tpd2(model,p,T,z,result,cache)
    comps, βi, volumes, idx_vapour = result
    phase_cache,ss_cache,newton_cache = cache
    pures,tpd_cache,lnϕv,lnϕl = phase_cache
    __cache,_,_,_,f1,f2,f3,_ = ss_cache
    di_cache,_,_ = __cache
    np = length(comps)
    nc = length(f1)
    lle = idx_vapour[] != 0
    zi,di,K = f1,f2,f3
    #calculate di,zi
    f1 .= 0
    f2 .= 0
    f3 .= 1
    for i in 1:np
        wi,βwi = comps[i],βi[i]
        f1 .+= βwi*wi
        dij = di_cache[i]
        lnϕwi,_ = lnϕ!(dij,model,p,T,wi,threaded = false,vol0 = volumes[i])
        dij .= lnϕwi .+ log.(wi)
        di .+= βwi*dij
    end
    lnϕ!(f3,model,p,T,z,threaded = false)
    dz = f3 .+ log.(z)
    dw = f3
    @show dz
    @show di
    for i in 1:np
        tpd_test = tpd(model,p,T,z,tpd_cache,break_first = false,lle = lle,di = di,strategy = :pure)
        @show tpd_test[1]
        @show tpd_test[2]
    end
    return K
end

function mm_di(model,p,T,result)
    comps, βi, volumes, idx_vapour = result
    di = similar(comps[1])
    for i in 1:1
        lnϕwi,_ = lnϕ!(di,model,p,T,comps[i],threaded = false,vol0 = volumes[i])
        di .= lnϕwi .+ log.(comps[i])
    end
    return di
end