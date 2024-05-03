struct MultiPhaseTPFlash{T} <: TPFlashMethod
    K0::Union{Vector{T},Nothing}
    ss_tol::Float64
    ss_iters::Int
    nacc::Int
    second_order::Bool
    max_phases::Int
    phase_iters::Int
end


function MultiPhaseTPFlash(;equilibrium = :vle,
    K0 = nothing,
    x0 = nothing,
    y0 = nothing,
    ss_tol = sqrt(eps(Float64)),
    ss_iters = 4,
    nacc = 5,
    second_order = false,
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

    return MultiPhaseTPFlash{T}(K0,ss_tol,ss_iters,nacc,second_order,max_phases,phase_iters)
end

function index_reduction(m::MultiPhaseTPFlash,idx::AbstractVector)
    K0,ss_tol,ss_iters,nacc,second_order,max_phases,phase_iters = m.K0,m.ss_tol,m.ss_iters,m.nacc,m.second_order,m.max_phases,m.phase_iters
    K0 !== nothing && (K0 = K0[idx])
    T = eltype(K0)
    return MultiPhaseTPFlash{T}(K0,ss_tol,ss_iters,nacc,second_order,max_phases,phase_iters)
end

function tpd_cache end

function tp_flash_multi_cache(model,p,T,z)
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
    v_cache = similar(f1)
    bi_cache = similar(f1)
    result_cache = comps_cache,bi_cache,v_cache
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

    #step 1: tpd: we check if there is any unstable phase.
    δn_add = _add_phases!(model,p,T,z,_result,phase_cache,options)
    if δn_add == 0
        return result
    end
    g0 = eos(model,volumes[1],T,z) + p*vz
    #step 2: main loop,iterate flashes until all phases are stable
    while !done
        iter += 1
        δn = δn_add - δn_remove
        if δn_add != 0 || δn_remove != 0#step 2.1: δn != 0. we add or remove phases and find new candidate ones.
            δn_add,δn_remove = 0,0
            n_phases += δn
            #sucessive substitution iteration
            cache = resize_cache!(cache,n_phases)
            result,ss_converged = tp_flash_multi_ss!(model,p,T,z,_result,ss_cache,options)
            δn_add = _add_phases!(model,p,T,z,_result,phase_cache,options)
            #@show δn_add
            if δn_add != 1
                δn_remove = _remove_phases!(model,p,T,z,_result,cache,options)
            end
            if δn_add == δn_remove == 0 && ss_converged
                @info "refinement reached! with RR"
                break
            end
        else #step 2.1: δn == 0. we refine our existing values.
            @info "refinement reached!"
            converged = true
            break
            #=
            comps, βi, volumes = tp_flash_multi_neq!(model,p,T,z,result,cache,options)
            #final stability check
            δn_remove = _remove_phases!(model,p,T,z,result,cache,options)
            δn_add = _add_phases!(model,p,T,z,result,cache,options)
            δn = δn_add - δn_remove
            δn == 0 && return result
            =#
        end
        done = iter > max_iter
        done = done || converged
    end
    gmix = _multiphase_gibbs(model,p,T,result)/(Rgas(model)*T)
    return comps, βi, volumes, gmix
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
    β = viewn(x0,np,np)
    β = view(β,1:np)
    β .= βi
    for i in 1:(np-1)
        lnKi = viewn(x0,np,i)
        xi = comps[i]
        lnKi .= log.(xi ./ xnp)
    end

    #options
    ss_iters = options.ss_iters
    ss_tol = options.ss_tol
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
        converged = dnorm(x0,x,1) < ss_tol
        converged && break
        x0 .= x
    end

    return result,converged
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
function multiphase_RR_β!(F, x, z, _result, ss_cache)
    _,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes, idx_vapour = _result
    result = comps, βi, volumes
    nc = length(z)
    np = length(βi)
    βF = viewn(F,np,np)

    #transforms lnK to K
    multiphase_lnK_K!(F_cache,x,np)
    if np == 2# && false
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
    #@show βsol
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

function multiphase_RR_lnK!(F, x, model, p, T, z, _result, ss_cache)
    _,_,_,F_cache,f1,f2,f3,dem_cache = ss_cache
    comps, βi, volumes, idx_vapour = _result
    result = comps, βi, volumes
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
    pures,tpd_cache = cache
    np = length(comps)
    nc = length(z)
    δn_add = 0
    max_phases = min(options.max_phases,nc)
    max_phases,np
    np + δn_add == max_phases && return 0 #we cannot add new phases here
    for i in 1:np
        #np + δn_add == max_phases && break
        w = comps[i]
        vw = volumes[i]
        diff_test = VT_diffusive_stability(model,vw,T,w)
        if !diff_test# && idx_vapour == 0
            δn_add += 1
            np += 1
            β1,x1,v1,β2,x2,v2 = split_phase_k(model,p,T,w,vw,pures,options.K0)
            β0 = β[i]
            β[i] = β0*β1
            comps[i] = x1
            volumes[i] = v1
            if pip(model,v2,T,x2) <= 1 && idx_vapour == 0
                _idx_vapour[] = np + 1
            end
            push!(comps,x2)
            push!(volumes,v2)
            push!(β,β0*β2)
        end
        is_lle = idx_vapour != i && idx_vapour != 0 #if we have a vapour phase, we only search for liquid-liquid splits
        tpd_i = tpd(model,p,T,w,tpd_cache,break_first = true, strategy = :pure,lle = is_lle)
        tpd_test = length(tpd_i) == 0
        if !tpd_test
            δn_add += 1
            np += 1
            tpd_comps,_,phases_z,phases_w = tpd_i
            #check our current phases and the trial ones
            y = tpd_comps[1]
            phase_w = phases_z[1]
            phase_y = phases_w[1]
            vy = volume(model,p,T,y,phase = phase_y)
            β0 = β[i]
            #phase not stable: generate a new one from tpd result
            β1,x1,v1,β2,x2,v2 = split_phase_tpd(model,p,T,w,y,phase_w,phase_y,vw)

            if is_vapour(phase_w) && idx_vapour == 0 #we identified the vapour phase
                _idx_vapour[] = i
            end
            β0 = β[i]
            β[i] = β0*β1
            comps[i] = x1
            volumes[i] = v1
            push!(comps,x2)
            push!(volumes,v2)
            push!(β,β0*β2)
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

function split_phase_tpd(model,p,T,z,w,phase_z,phase_w,vz = volume(model,p,T,z,phase = phase_z),vw = volume(model,p,T,w,phase = phase_w))
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

    return βi,x3,v3,(1-βi),w,vw
end

function split_phase_k(model,p,T,z,vz,pures = split_model(model),K = nothing)
    #split phase via K values.
    if isnothing(K)
        Ki = suggest_K(model,p,T,z,pures)
    else
        Ki = K
    end
    βi,singlephase,_ = rachfordrice_β0(Ki,z)
    wx = similar(z)
    wy = similar(z)
    rr_flash_liquid!(wx,Ki,z,βi)
    wy .= wx .* Ki
    wx ./= sum(wx)
    wy ./= sum(wy)
    vx = volume(model,p,T,wx)
    vy = volume(model,p,T,wy)
    return 1-βi,wx,vx,βi,wy,vy
end
