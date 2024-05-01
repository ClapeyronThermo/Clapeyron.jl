function tpd_cache end

function tp_flash_michelsen_multi_cache(model,p,T,z)
    pure = split_model(model)
    _tpd_cache = tpd_cache(model,p,T,z)
    phase_cache = (pure,_tpd_cache)
    F_cache = zeros(Base.promote_eltype(model,p,T,z),length(z))
    F_cache2 = similar(F_cache)
    F_cache3 = similar(F_cache)
    f1,f2,f3 = similar(F_cache),similar(F_cache),similar(F_cache)
    F3,F4,F5,ΔF1,ΔF2,xdem,Fdem = similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache),similar(F_cache)
    dem_cache = (F3,F4,F5,ΔF1,ΔF2,xdem,Fdem)
    comps_cache = fill(similar(f1),1)
    ss_cache = comps_cache,F_cache,F_cache2,F_cache3,f1,f2,f3
    return phase_cache,ss_cache,nothing
end

function resize_cache!(cache,np)
    phase_cache,ss_cache,newton_cache = cache
    nc = length(phase_cache[1])
    comps_cache,F_cache,F_cache2,F_cache3,f1,f2,f3,dem_cache = ss_cache
    np_old = length(comps_cache)
    resize!(comps_cache,np)
    for i in (np_old+1):np
        comps_cache[i] = similar(f1)
    end
    for fi in (F_cache,F_cache2,F_cache3)
        resize!(fi,nc*np)
    end

    for fi in dem_cache
        resize!(fi,nc*np)
    end
    return cache
end

function tp_flash_michelsen_multi(model,p,T,z,options = nothing)
    n_phases = 1
    max_iter = 20
    #step 1: tpd: we check if there is any unstable phase.
    comps,βi,_,_ = tpd(model,p,T,z,strategy = :pure,break_first = true)
    #step 0: no phases found.
    if !iszero(length(comps))
        resize!(comps,0)
        resize!(βi,0)
    end

    push!(comps,z ./ sum(z))
    push!(βi,1.0)
    volumes = [volume(model,p,T,z)]

    iszero(length(comps)) && return comps,βi,[v]
    result = (comps, βi, volumes)
    δn_add = 1
    δn_remove = 0

    cache = tp_flash_michelsen_multi_cache(model,p,T,z)
    phase_cache,ss_cache,newton_cache = cache
    options = nothing
    done = false
    converged = false
    iter = 0

    δn_add =_add_phases!(model,p,T,z,result,phase_cache,options)
    if δn_add == 0
        return result
    end
    #step 2: main loop,iterate flashes until all phases are stable
    while !done
        iter += 1
        δn = δn_add - δn_remove
        @show δn_remove,δn_add
        if δn_add != 0 || δn_remove != 0#step 2.1: δn != 0. we add or remove phases and find new candidate ones.

            n_phases += δn

            #TODO: add a limit on how many phases exist

            #sucessive substitution iteration
            cache = resize_cache!(cache,n_phases)
            comps, βi, volumes = tp_flash_michelsen_n_ss!(model,p,T,z,result,ss_cache,options)
            return result
            #δn_remove = _remove_phases!(model,p,T,z,result,phase_cache,options)
            δn_add = _add_phases!(model,p,T,z,result,phase_cache,options)
            if δn_add == 1
                δn_remove = 0
            else
                δn_remove = _remove_phases!(model,p,T,z,result,cache,options)
            end
            @show δn_remove,δn_add
        else #step 2.1: δn == 0. we refine our existing values.
            @info "refinement reached!"
            converged = true
            #=
            comps, βi, volumes = tp_flash_michelsen_n_neq!(model,p,T,z,result,cache,options)
            #we check the diffusive stability of each phase.
            #final stability check
            δn_remove = _remove_phases!(model,p,T,z,result,cache,options)
            δn_add = _add_phases!(model,p,T,z,result,cache,options)
            δn = δn_add - δn_remove
            δn == 0 && return result
            =#
            return result
        end
        done = iter > max_iter
        @show done
        done = done || converged
        @show converged
    end
    
    return result
end

function tp_flash_michelsen_n_ss!(model,p,T,z,result,ss_cache,options)
    comps_cache,x0,x,xx,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes = result
    nc = length(z)
    np = length(βi)
    #set the initial values to 0
    x0 .= 0
    x .= 0
    xx .= 0
    βi_cache = view(f1,1:np)
    volumes_cache = view(f2,1:np)
    result_cache = comps_cache,βi_cache,volumes_cache
    #fill F with data
    xnp = comps[np]
    β = viewn(x0,np,np)
    β = view(β,1:np)
    β .= βi
    for i in 1:(np-1)
        lnKi = viewn(x0,np,i)
        xi = comps[i]
        lnKi .= log.(xi ./ xnp)
    end
    
    #options
    max_iters = 10*n_phases
    itacc = 0
    nacc = 5
    ss_tol = 1e-8
    
    F3,F4,F5,ΔF1,ΔF2,xdem,Fdem = dem_cache
    for i in 1:max_iters
        itacc += 1
        fixpoint_multiphase!(x, x0, model, p, T, z, result, ss_cache)
        #store values in cache
        βi_cache .= βi
        volumes_cache .= volumes
        for i in 1:np
            xi_cache,xi = comps_cache[i],comps[i]
            xi_cache .= xi
        end
        
        #acceleration
        if itacc == (nacc - 2)
            F3 .= F
        elseif itacc == (nacc - 1)
            F4 .= F
        elseif itacc == nacc
            itacc = 0
            F5 .= F
            # acceleration using DEM (1 eigenvalues)
            xdem = dem!(xdem, F5, F4, F3,(ΔF1,ΔF2)) 
            fixpoint_multiphase!(Fdem, xdem, model, p, T, z, result_cache, ss_cache)
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

        dnorm(x0,x,1) < K_tol && break
        x0 .= x
    end

    return result
end

function fixpoint_multiphase!(F, x, model, p, T, z, result, ss_cache)
    comps_cache,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes = result
    #given constant K, calculate β
    multiphase_RR_β!(F, x, z, result, ss_cache)
    #given constant β, calculate K
    multiphase_RR_lnK!(F, x, model, p, T, z, result, ss_cache)
    return F
end

function RR_t!(t,x,β,np,nc)
    t .= 1
    for l in 1:(np - 1)
        Kl = viewn(x,np,l)
        βl = β[l]
        fi = zero(βl)
        for i in 1:nc
            #K at phase i
            t[i] += βl*expm1(Kl[i])
        end
    end
    return t
end

#constant K, calculate β
function multiphase_RR_β!(F, x, z, result, ss_cache)
    comps_cache,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes = result
    nc = length(z)
    np = length(βi)
    βF = viewn(F,np,np)
    
    #transforms lnK to K
    multiphase_lnK_K!(F_cache,x,np)
    if np == 2 && false
        K = viewn(F_cache,np,1)
        βsol = rachfordrice(K,z)
        βF[1] = βsol #two-phase solution does not require solving a neq problem
        βF[2] = 1 - βsol
        βi[1] = βsol
        βi[2] = 1 - βsol
        return F
    end

    function f(β)
        #okuno et al. (2010): objetive function
        t = similar(β,nc)
        t = RR_t!(t,x,β,np,nc) #we use lnK here.
        res = zero(eltype(β))
        for i in 1:nc
            res -= z[i]*log(abs(t[i]))
        end
        return res
    end
    β0 = viewn(x,np,np)
    β0 = copy(view(β0,1:np-1))
    result = Solvers.optimize(f,β0,LineSearch(Newton(linsolve = static_linsolve)))
    βsol = Solvers.x_sol(result)
    ∑β = sum(βsol)
    for i in 1:(np-1)
        βF[i] = βsol[i]
        βi[i] = βsol[i]
    end
    @show βsol
    #update last phase fraction
    βF[np] = 1 - ∑β
    βi[np] = 1 - ∑β
    return F
end

function multiphase_lnK_K!(F,x,np)
    for l in 1:(np - 1)
        lnK = viewn(x,np,l)
        K = viewn(F,np,l)
        K .= exp.(lnK)
    end
    return F
end

function multiphase_K_lnK!(F,x,np)
    for l in 1:(np - 1)
        K = viewn(x,np,l)
        lnK = viewn(F,np,l)
        lnK .= log.(K)
    end
    return F
end

function multiphase_RR_lnK!(F, x, model, p, T, z, result, ss_cache)
    comps_cache,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes = result
    nc = length(z)
    np = length(βi)
    multiphase_lnK_K!(F_cache,x,np) #F_cache now contains Kij
    β = viewn(F,np,np)
    
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
        Ki = viewn(F_cache,np,i)
        xi = comps[i]
        xi .= Ki .* z ./ t
        xi ./= sum(xi)
    end
    
    #update lnK
    for i in 1:np - 1
        lnKi = viewn(F,np,i)
        xi = comps[i]
        lnϕi,vi = lnϕ!(f3, model, p, T, xi, vol0 = volumes[i])
        isnan(vi) && (lnϕi,vi = lnϕ!(f3, model, p, T, xi)) #restart
        volumes[i] = vi
        lnKi .=  lnϕnp .- lnϕi
    end
    @show viewn(F,np,1)
    return F
end
#perform tpd on each phase, if not stable, generate a new one
function _add_phases!(model,p,T,z,result,cache,options)
    comps, β, volumes = result
    pures,tpd_cache = cache
    n = length(comps)
    np = n
    nc = length(z)
    δn_add = 0
    for i in 1:n
        np == nc && continue
        w = comps[i]
        tpd_test = tpd(model,p,T,w,tpd_cache,break_first = true, strategy = :pure)
        if length(tpd_test[1]) > 0
            δn_add += 1
            np += 1
            #phase not stable: generate a new one:
            Ki = suggest_K(model,p,T,w,pures)

            #TODO: RR refine here ?
            βi,singlephase,_ = rachfordrice_β0(Ki,w)
            wx = similar(w)
            wy = similar(w)
            rr_flash_liquid!(wx,Ki,z,βi)
            wy .= wx .* Ki
            wx ./= sum(wx)
            wy ./= sum(wy)
            vx = volume(model,p,T,wx)
            vy = volume(model,p,T,wy)
            β0 = β[i]
            β[i] = β0*βi
            comps[i] = wy
            volumes[i] = vy
            push!(comps,wx)
            push!(volumes,vx)
            push!(β,β0*(1 - βi))
        end
    end
    return δn_add
end

#find the phase with the minimum βi. if mixing that phase with any other phase generates
#a more stable phase, remove it
function _remove_phases!(model,p,T,z,result,cache,options)
    comps, β, volumes = result
    n = length(comps)
    δn_remove = 0
    comps, β, volumes = result
    pures,tpd_cache = cache
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
        gi = VT_gibbs_free_energy(model,vi,T,wi)
        gmix = VT_gibbs_free_energy(model,vmix,T,wmix)
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

function _multiphase_gibbs(model,p,T,result)
    comps, β, volumes = result
    res = zero(Base,promote_eltype(model,p,T,comps[1]))
    np = length(comps)
    for i in 1:np
        res += β[i] * VT_gibbs_free_energy(model,volumes[i],T,comps[i])
    end
    return res
end

