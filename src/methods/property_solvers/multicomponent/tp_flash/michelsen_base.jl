function michelsen_optimization_of!(g,H,model,p,T,z,caches,ny_var,gz)
    second_order = !isnothing(H)
    nx,ny,vcache,lnϕ_cache,in_eq,phases = caches
    in_equilibria,non_inx,non_iny = in_eq
    phasex,phasey = phases
    volx,voly = vcache[]
    update_nxy!(nx,ny,ny_var,z,non_inx,non_iny) #updates nx, ny with ny_var vector
    nxsum = sum(nx)
    nysum = sum(ny)
    x = nx
    y = ny
    nx .= nx ./ nxsum
    ny .= ny ./ nysum
    f = zero(eltype(ny_var))
    f -= gz
    !isnothing(g) && (g .= 0)
    if second_order
        lnϕx, ∂lnϕ∂nx, volx = modified_∂lnϕ∂n(model, p, T, x, lnϕ_cache; phase=phasex, vol0=volx)
        ∂x,∂2x = lnϕx,∂lnϕ∂nx
        ∂2x .-= 1
        ∂2x ./= nxsum
        for i in axes(∂2x,1)
            ∂2x[i,i] += 1/(nxsum*x[i])
            ∂x[i] += log(x[i])
            non_inx[i] && (∂x[i] = 0)
        end

        H .= @view ∂2x[in_equilibria, in_equilibria]
        !isnothing(g) && (g .= @view ∂x[in_equilibria])
        f += nxsum*dot(∂x,x)

        lnϕy, ∂lnϕ∂ny, voly = modified_∂lnϕ∂n(model, p, T, y, lnϕ_cache; phase=phasey, vol0=voly)
        ∂y,∂2y = lnϕy,∂lnϕ∂ny
        ∂2y .-= 1
        ∂2y ./= nysum
        for i in axes(∂2y,1)
            ∂2y[i,i] += 1/(nysum*y[i])
            ∂y[i] += log(y[i])
            non_iny[i] && (∂y[i] = 0)
        end

        H .+= @view ∂2y[in_equilibria, in_equilibria]
        !isnothing(g) && (g .-= @view ∂y[in_equilibria])
        !isnothing(g) && (g .*= -1)

        f += nysum*dot(∂y,y)
    else
        ∂x,volx = modified_lnϕ(model, p, T, x,lnϕ_cache; phase=phasex, vol0=volx)
        for i in axes(∂x,1)
            ∂x[i] += log(x[i])
            non_inx[i] && (∂x[i] = 0)
        end
        !isnothing(g) && (g .= @view ∂x[in_equilibria])
        f += nxsum*dot(∂x,x)
        ∂y,voly = modified_lnϕ(model, p, T, y, lnϕ_cache; phase=phasey, vol0=voly)
        for i in axes(∂y,1)
            ∂y[i] += log(y[i])
            non_iny[i] && (∂y[i] = 0)
        end
        !isnothing(g) && (g .-= @view ∂y[in_equilibria])
        !isnothing(g) && (g .*= -1)
        f += nysum*dot(∂y,y)
    end
    vcache[] = (volx,voly)
    return f
end

function michelsen_gibbs_feed(model,p,T,z,caches)
    nx,ny,vcache,lnϕ_cache,in_eq,phases = caches
    in_equilibria,non_inx,non_iny = in_eq
    if all(in_equilibria) && has_a_res(model)
        ∂z,volz = modified_lnϕ(model, p, T, z, lnϕ_cache)
        ∂z .+= log.(z)
        gz = dot(@view(z[in_equilibria]),@view(∂z[in_equilibria]))
    else
        gz = zero(Base.promote_eltype(model,p,T,z))
    end
    return gz
end

function michelsen_optimization_obj(model,p,T,z,caches)
    gz = michelsen_gibbs_feed(model,p,T,z,caches)

    function objective_ip(x)
        return michelsen_optimization_of!(nothing,nothing,model,p,T,z,caches,x,gz)
    end

    function gradient_ip(∇f, x)
        michelsen_optimization_of!(∇f,nothing,model,p,T,z,caches,x,gz)
        return ∇f
    end

    function objective_gradient_ip(∇f, x)
        fx = michelsen_optimization_of!(∇f,nothing,model,p,T,z,caches,x,gz)
        return fx,∇f
    end
    function hessian_ip(∇²f, x)
        michelsen_optimization_of!(nothing,∇²f,model,p,T,z,caches,x,gz)
        return ∇²f
    end
    function objective_gradient_hessian_ip(∇f, ∇²f, x)
        fx = michelsen_optimization_of!(∇f,∇²f,model,p,T,z,caches,x,gz)
        return fx, ∇f, ∇²f
    end

    scalarobj_ip = NLSolvers.ScalarObjective(f=objective_ip,
                                g=gradient_ip,
                                fg=objective_gradient_ip,
                                fgh=objective_gradient_hessian_ip,
                                h=hessian_ip)
    return scalarobj_ip
    #optprob_ip = NLSolvers.OptimizationProblem(scalarobj_ip; inplace=true)
end


#updates lnK, returns lnK,volx,voly, gibbs if β !== nothing
function update_K!(lnK,model,p,T,x,y,z,β,vols,phases,non_inw,dlnϕ_cache = nothing)
    volx,voly = vols
    phasex,phasey = phases
    non_inx,non_iny = non_inw
    lnϕx, volx = modified_lnϕ(model, p, T, x, dlnϕ_cache; phase = phasex, vol0=volx)

    lnK .= lnϕx
    gibbs = zero(eltype(lnK))
    if β !== nothing
        for i in eachindex(y)
            if !non_inx[i] || isinf(lnK[i])
                 gibbs += (1-β)*x[i]*(log(x[i]) + lnϕx[i])
            end
        end
    else
        gibbs = gibbs/gibbs
    end

    lnϕy, voly = modified_lnϕ(model, p, T, y, dlnϕ_cache; phase = phasey, vol0=voly)
    lnK .-= lnϕy
    if β !== nothing
        for i in eachindex(y)
            if !non_iny[i] || iszero(exp(lnK[i]))
                gibbs += β*y[i]*(log(y[i]) + lnϕy[i])
            end
        end
    else
        gibbs = gibbs/gibbs
    end

    return lnK,volx,voly,gibbs
end

#updates x,y after a sucessful rachford rice procedure
function update_rr!(K,β,z,x,y,
                    non_inx=nothing,
                    non_iny=nothing;normalize::Bool = true)

    RR = RachfordRiceResidual(K,z,non_inx,non_iny)
    result = RachfordRiceResult(β)
    return update_rr!(x,y,RR,result;normalize)
end

function update_rr!(x,y,K::RachfordRiceResidual,res::RachfordRiceResult;normalize::Bool = true)
    z = K.z
    _1 = one(eltype(x))
    _0 = zero(_1)
    βv,βl = res.V,res.L
    if normalize
        n = sum(z)
    else
        n = one(eltype(z))
    end

    for i in eachindex(x)
        Ki,zi = K[i],z[i]
        wi = zi/n
        if isinf(Ki)
            x[i] = _0
            y[i] = wi/βv
        else
            wi = zi/n
            detKi = RachfordRice.inv_fma(Ki,βv,βl)
            xi = detKi*wi
            x[i] = xi
            y[i] = xi*Ki
        end
    end
    return x,y
end

function update_nxy!(nx,ny,ny_var,z,non_inx,non_iny)
    ii = 0
    for i in eachindex(z)
        if non_inx[i]
            ny[i] = z[i]
            nx[i] = 0.0
        elseif non_iny[i]
            ny[i] = 0.0
            nx[i] = z[i]
        else
            ii += 1
            nyi = ny_var[ii]
            ny[i] = nyi
            nx[i] = z[i] - nyi
        end
    end
    return nx,ny
end

function tp_flash_K0(model,p,T,z)
    K = zeros(Base.promote_eltype(model,p,T,z),length(model))
    tp_flash_K0!(K,model,p,T,z,nothing)
    return K
end

function tp_flash_K0!(K,model,p,T,z,cache)
    K_calculated = tp_flash_fast_K0!(K,model,p,T,z)

    if K_calculated
        Kmin,Kmax = extrema(K)
        if Kmin >= 1 || Kmax <= 1
            K_calculated = false
        end
    end

    if !K_calculated
        suggest_K!(K,model,p,T,z,cache)
    end
end

function tp_flash_fast_K0!(K,model,p,T,z)
    return false
end

function pt_flash_x0(model,p,T,n,method = GeneralizedXYFlash(),non_inx = FillArrays.Fill(false,length(model)),non_iny = FillArrays.Fill(false,length(model)))
    ∑n = sum(n)
    z = n/∑n
    TT = Base.promote_eltype(model,p,T,n)
    if is_vle(method)
        phasex,phasey = :liquid,:vapour
    elseif is_lle(method)
        phasex,phasey = :liquid,:liquid
    else
        phasex,phasey = :unknown,:unknown
    end
    non_inw = (non_inx,non_iny)
    _1,_0 = one(TT),zero(TT)
    x,y = similar(z,TT),similar(z,TT)
    x .= z
    y .= z
    K,lnK = similar(z,TT),similar(z,TT)
    verbose = get_verbosity(method)
    RR = RachfordRiceResidual(K,z;non_inx,non_iny)
    if !isnothing(method.K0)
        K .= _1 * method.K0
        lnK .= log.(K)
        volx = zero(_1)
        voly = zero(_1)
        verbose && @info "K0 already provided"
    elseif hasfield(typeof(method),:x0) && hasfield(typeof(method),:y0) && !isnothing(method.x0) && !isnothing(method.y0)
        x .= method.x0 ./ sum(method.x0)
        y .= method.y0 ./ sum(method.y0)
        lnK .= log.(x ./ y)
        volx = zero(_1)
        voly = zero(_1)
        if method.v0 === nothing
            lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,z,nothing,(nothing,nothing),(phasex,phasey),non_inw)
        else
            vl0,vv0 = method.v0
            lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,z,nothing,(vl0,vv0),(phasex,phasey),non_inw)
        end
        K .= exp.(lnK)
        verbose && @info "x0,y0 provided, calculating K0 via Clapeyron.update_K!"
    elseif is_vle(method) || is_unknown(method)
        # VLE correlation for K
        verbose && @info "K0 calculated via pure VLE correlation"
        tp_flash_K0!(K,model,p,T,z,nothing)
        phasex = :liquid
        phasey = :vapour

        #if we can't predict K, we use lle
        if is_unknown(method)

            if !rr_feasible(RR)
                verbose && @info "VLE correlation failed, trying LLE initial point."
                K_lle = K0_lle_init(model,p,T,z;reduced = true)
                if any(!isone,K_lle) #only use LLE result if actually exists
                    K .= K_lle
                    phasex = :liquid
                    phasey = :liquid
                end
                lnK .= log.(K)
            else
                phasex,phasey = :liquid,:vapour
            end
        end
        lnK .= log.(K)
        volx = zero(_1)
        voly = zero(_1)
    else
        verbose && @info "K0 calculated via LLE initial point (tpd)"
        K .= K0_lle_init(model,p,T,z;reduced = true)
        lnK .= log.(K)
        phasey = :liquid
        volx = zero(_1)
        voly = zero(_1)
    end

    verbose && @info "K0 = $K"
    # Initial guess for phase split
    status = rr_status(RR,K_tol = 1e-10)

#=
    TREND bubble/dew initialization
    Maybe initial K values overshoot the actual phase split.
    if the initial K values generate a single phase result, but we can split the K into two compositions (Kmin < 1 or Kmax > 1)
    then we start at the bubble (or dew conditions)
    =#

    if status == RRLiquid || status == RRVapour
        bubble = status == RRLiquid
        if !rr_feasible(RR)
            β = rr_feasibility_projection!(RR,z,bubble;verbose = verbose)
            x,y = update_rr!(K,β,z,x,y,non_inx,non_iny)
        else
            β = bubble ? _0 : _1
        end
        status = RREq
        βsgn = ifelse(bubble,1,-1)
        β += βsgn*eps(eltype(β))
        β_result = RachfordRiceResult(β)
    elseif status == RREq
        β_result = rachfordrice(RR; K_tol, verbose)
        β = β_result.V
    else
        β = _0/_0
        β_result = RachfordRiceResult(β)
    end

    verbose && @info "initial vapour fraction = $β"
    #verbose && status != RREq && @info "initial point is single-phase (does not satisfy Rachford-Rice constraints). Exiting early"

    update_rr!(x,y,RR,β_result;normalize = true)
    βv = ∑n*β_result.V
    βl = ∑n*β_result.L
    if hasfield(typeof(method),:v0) && !isnothing(method.v0) && iszero(volx) && iszero(voly)
        vl0,vv0 = method.v0
        volx,voly = _1*vl0,_1*vv0
    end

    iszero(volx) && (volx = volume(model,p,T,x,phase = phasex))
    iszero(voly) && (voly = volume(model,p,T,y,phase = phasey))
    lle = is_liquid(phasex) && is_liquid(phasey)
    vapour_idx = lle ? -1 : 2
    data = FlashData(p,T,zero(T),vapour_idx)
    #has_a_res(model) && is_liquid(VT_identify_phase(model,voly,T,y)) && (voly = Rgas(model)*T/p)
    r = FlashResult(SA[x,y],SA[βl,βv],SA[volx,voly],data)
    return r
end

#=
function dgibbs_obj!(model::EoSModel, p, T, z, phasex, phasey,
    nx, ny, vcache, ny_var = nothing, in_equilibria = FillArrays.Fill(true,length(z)), non_inx = in_equilibria, non_iny = in_equilibria;
    F=nothing, G=nothing, H=nothing)

    # Objetive Function to minimize the Gibbs energy
    # It computes the Gibbs energy, its gradient and its hessian
    iv = 0
    for i in eachindex(z)
        if in_equilibria[i]
            iv += 1
            nyi = ny_var[iv]
            ny[i] = nyi
            nx[i] = z[i] - nyi
        end
    end    # nx = z .- ny

    nxsum = sum(nx)
    nysum = sum(ny)
    x = nx ./ nxsum
    y = ny ./ nysum

    # Volumes are set from local cache to reuse their values for following
    # Iterations
    volx,voly = vcache[]
    all_equilibria = all(in_equilibria)
    if H !== nothing
        # Computing Gibbs energy Hessian
        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y; phase=phasey, vol0=voly)

        if !all_equilibria
            ∂ϕx = ∂lnϕ∂nx[in_equilibria, in_equilibria]
            ∂ϕy = ∂lnϕ∂ny[in_equilibria, in_equilibria]
        else
            #skip a copy if possible
            ∂ϕx,∂ϕy = ∂lnϕ∂nx,∂lnϕ∂ny
        end
            ∂ϕx .-= 1
            ∂ϕy .-= 1
            ∂ϕx ./= nxsum
            ∂ϕy ./= nysum
        for (i,idiag) in pairs(diagind(∂ϕy))
            ∂ϕx[idiag] += 1/nx[i]
            ∂ϕy[idiag] += 1/ny[i]
        end

        #∂ϕx = eye./nx .- 1/nxsum .+ ∂lnϕ∂nx/nxsum
        #∂ϕy = eye./ny .- 1/nysum .+ ∂lnϕ∂ny/nysum
        H .= ∂ϕx .+ ∂ϕy
    else
        lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
    end
    #volumes are stored in the local cache
    vcache[] = (volx,voly)

    ϕx = log.(x) .+ lnϕx #log(xi*ϕxi)
    ϕy = log.(y) .+ lnϕy #log(yi*ϕyi)

    # to avoid NaN in Gibbs energy
    for i in eachindex(z)
        non_iny[i] && (ϕy[i] = 0.)
        non_inx[i] && (ϕx[i] = 0.)
    end

    if G !== nothing
        # Computing Gibbs energy gradient
        i0 = 0
        for i in eachindex(in_equilibria)
            if in_equilibria[i]
                i0 += 1
                G[i0] = ϕy[i] - ϕx[i]
            end
        end
    end

    if F !== nothing
        # Computing Gibbs energy
        FO = dot(ny,ϕy) + dot(nx,ϕx)
        return FO
    end
end
=#

