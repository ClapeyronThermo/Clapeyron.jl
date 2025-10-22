function rachfordrice(K, z; β0=nothing,K_tol = 4*eps(eltype(K)), non_inx=FillArrays.Fill(false,length(z)), non_iny=FillArrays.Fill(false,length(z)))
    # Function to solve Rachdord-Rice mass balance
    β,status,limits = rachfordrice_β0(K,z,β0,non_inx,non_iny;K_tol = K_tol)
    if length(z) <= 3 && all(Base.Fix2(>,0),z) && all(!,non_inx) && all(!,non_iny) && status == RREq
        βx = rr_vle_vapor_fraction_exact(K,z)
        return clamp(βx,zero(β),one(β))
    end

    if status == RREq
        βx = rr_flash_refine(K, z, β, non_inx, non_iny, limits) # bracketed Halley when possible
        return clamp(βx,zero(β),one(β))
    elseif status == RRLiquid
        return zero(β)   # or eps(eltype(β))
    elseif status == RRVapour
        return one(β)   # or 1 - eps(eltype(β))
    else
        return zero(β)/zero(β)
    end
end

function K_extrema(K::AbstractVector{T},non_inx,non_iny) where T
    Kmax = T(-Inf)
    Kmin = T(Inf)
    for i in 1:length(K)
        if non_inx[i]
            Ki = T(Inf)
        elseif non_iny[i]
            Ki = T(0)
        else
            Ki = K[i]
        end
        Kmax = max(Ki,Kmax)
        Kmin = min(Ki,Kmin)
    end
    return Kmin,Kmax
end

function rr_margin_check(K,z,non_inx = FillArrays.Fill(false,length(K)),non_iny = FillArrays.Fill(false,length(K));K_tol = sqrt(eps(eltype(K))),verbose = false)
    Ktype = eltype(K)
    _0 = zero(Ktype)
    _1 = one(Ktype)
    F0 = Clapeyron.rr_flash_eval(K, z, _0, non_inx, non_iny) # F(0)
    F1 = Clapeyron.rr_flash_eval(K, z, _1, non_inx, non_iny) # F(1)
    cond0 = (abs(F0) <= K_tol) && (F1 < 0) # bubble candidate
    cond1 = (abs(F1) <= K_tol) && (F0 > 0) # dew candidate
    
    if verbose
        δβ = sqrt(eps(Ktype))
        Fp0_num = (Clapeyron.rr_flash_eval(K, z, δβ, non_inx, non_iny) - F0) / δβ
        @info """ checking boundary conditions: 
        F(0)             = $(F0)
        F(1)             = $(F1) 
        F'(0)            ≈ $(Fp0_num) 
        gate             = $(cbrt(K_tol))
        bubble condition = $(cond0) 
        dew condition    = $(cond1)
        """
    end

    if cond0 && !cond1
        β = _0
        status = RRLiquid
        verbose && @info "boundary projection applied: β→0 (bubble), status=RRLiquid"
    elseif !cond0 && cond1
        β = _1
        status = RRVapour
        verbose && @info "boundary projection applied: β→1 (dew), status=RRVapour"
    elseif cond0 && cond1
        # choose the boundary with smaller linearized step δβ = |F|/|F'|
        Fp0,Fp1 = _0,_0
        @inbounds for i in eachindex(K,z)
            Δ   = K[i]-_1
            zi  = z[i]
            Fp0 -= zi*(Δ*Δ)
            Fp1 -= zi*(Δ*Δ)/(K[i]*K[i])
        end
        δβ0 = abs(F0)/(abs(Fp0)+eps(Ktype))
        δβ1 = abs(F1)/(abs(Fp1)+eps(Ktype))
        if δβ0 <= δβ1
            β = _0
            status = RRLiquid
            verbose && @info "boundary projection via linearized distance: β→0 (δβ0=$(δβ0) ≤ δβ1=$(δβ1))"
        else
            β = _1
            status = RRVapour
            verbose && @info "boundary projection via linearized distance: β→1 (δβ1=$(δβ1) < δβ0=$(δβ0))"
        end
    else
        verbose && @info "boundary projection not conclusive."
        status = RREq
        β = _0/_0
    end
    return status,β
end

function michelsen_optimization_of!(g,H,model,p,T,z,caches,ny_var,gz)
    second_order = !isnothing(H)
    nx,ny,vcache,lnϕ_cache,in_eq,phases = caches
    in_equilibria,non_inx,non_iny = in_eq
    phasex,phasey = phases
    volx,voly = vcache[]
    update_nxy!(nx,ny,ny_var,z,non_inx,non_iny) #updates nx, ny with ny_var vector
    nxsum = sum(nx)
    nysum = sum(ny)
    x = nx ./ nxsum
    y = ny ./ nysum
    f = zero(eltype(ny_var))
    f -= gz
    !isnothing(g) && (g .= 0) 
    if second_order
        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x, lnϕ_cache; phase=phasex, vol0=volx)
        ∂x,∂2x = lnϕx,∂lnϕ∂nx
        ∂2x .-= 1
        ∂2x ./= sum(nx)
        for i in 1:size(∂2x,1)
            ∂2x[i,i] += 1/nx[i]
            ∂x[i] += log(x[i])
            non_inx[i] && (∂x[i] = 0)
        end

        H .= @view ∂2x[in_equilibria, in_equilibria]
        !isnothing(g) && (g .= @view ∂x[in_equilibria])
        f += dot(∂x,nx)

        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y, lnϕ_cache; phase=phasey, vol0=voly)
        ∂y,∂2y = lnϕy,∂lnϕ∂ny
        ∂2y .-= 1
        ∂2y ./= sum(ny)
        for i in 1:size(∂2y,1)
            ∂2y[i,i] += 1/ny[i]
            ∂y[i] += log(y[i])
            non_iny[i] && (∂y[i] = 0)
        end

        H .+= @view ∂2y[in_equilibria, in_equilibria]
        !isnothing(g) && (g .-= @view ∂y[in_equilibria])
        !isnothing(g) && (g .*= -1)

        f += dot(∂y,ny)
    else
        ∂x,volx = lnϕ(model, p, T, x,lnϕ_cache; phase=phasex, vol0=volx)
        for i in 1:size(∂x,1)
            ∂x[i] += log(x[i])
            non_inx[i] && (∂x[i] = 0)
        end
        !isnothing(g) && (g .= @view ∂x[in_equilibria])
        f += dot(∂x,nx)
        ∂y,voly = lnϕ(model, p, T, y,lnϕ_cache; phase=phasey, vol0=voly)
        for i in 1:size(∂y,1)
            ∂y[i] += log(y[i])
            non_iny[i] && (∂y[i] = 0)
        end
        !isnothing(g) && (g .-= @view ∂y[in_equilibria])
        !isnothing(g) && (g .*= -1)
        f += dot(∂y,ny)
    end
    vcache[] = (volx,voly)
    return f
end

function michelsen_gibbs_feed(model,p,T,z,caches)
    nx,ny,vcache,lnϕ_cache,in_eq,phases = caches
    in_equilibria,non_inx,non_iny = in_eq
    if all(in_equilibria)
        ∂z,volz = lnϕ(model, p, T, z, lnϕ_cache)
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
        fx = michelsen_optimization_of!(nothing,nothing,model,p,T,z,caches,x,gz)
    end

    function gradient_ip(∇f, x)
        fx = michelsen_optimization_of!(∇f,nothing,model,p,T,z,caches,x,gz)
        return ∇f
    end

    function objective_gradient_ip(∇f, x) 
        fx = michelsen_optimization_of!(∇f,nothing,model,p,T,z,caches,x,gz)
        return fx,∇f
    end
    function hessian_ip(∇²f, x)
        fx = michelsen_optimization_of!(nothing,∇²f,model,p,T,z,caches,x,gz)
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
    #optprob_ip = NLSolvers.OptimizationProblem(scalarobj_ip; inplace=true)
end


#updates lnK, returns lnK,volx,voly, gibbs if β != nothing
function update_K!(lnK,model,p,T,x,y,z,β,vols,phases,non_inw,dlnϕ_cache = nothing)
    volx,voly = vols
    phasex,phasey = phases
    non_inx,non_iny = non_inw
    lnϕx, volx = lnϕ(model, p, T, x, dlnϕ_cache; phase = phasex, vol0=volx)
    if isnan(volx)
        lnϕx, volx = lnϕ(model, p, T, x, dlnϕ_cache, phase = phasex)
    end
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

    lnϕy, voly = lnϕ(model, p, T, y, dlnϕ_cache; phase = phasey, vol0=voly)
    if isnan(voly)
        lnϕy, voly = lnϕ(model, p, T, y, dlnϕ_cache, phase = phasey)
    end
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
                    non_inx=FillArrays.Fill(false,length(z)),
                    non_iny=FillArrays.Fill(false,length(z)))

    x = rr_flash_liquid!(x,K,z,β)
    y .= x .* K
    for i in eachindex(z)
        # modification for non-in-y components Ki -> 0
        if non_iny[i] || iszero(K[i])
            x[i] = z[i] / (1. - β)
            y[i] = 0.
        end
        # modification for non-in-x components Ki -> ∞
        if non_inx[i] || isinf(K[i])
            x[i] = 0.
            y[i] = z[i] / β
        end
    end
    x ./= sum(x)
    y ./= sum(y)
    return x,y
end

function update_nxy!(nx,ny,ny_var,z,non_inx,non_iny)
    iv = 0mar
    for i in eachindex(z)            
        if non_inx[i]
            ny[i] = z[i]
            nx[i] = 0.0
        elseif non_iny[i]
            ny[i] = 0.0
            nx[i] = z[i]
        else
            iv += 1
            nyi = ny_var[iv]
            ny[i] = nyi
            nx[i] = z[i] - nyi
        end
    end 
    return nx,ny
end

function tp_flash_K0(model,p,T,z)
    K = zeros(Base.promote_eltype(model,p,T,z),length(model))
    tp_flash_K0!(K,model,p,T,z)
    return K
end

function tp_flash_K0!(K,model,p,T,z)
    K_calculated = tp_flash_fast_K0!(K,model,p,T,z)

    if K_calculated
        Kmin,Kmax = extrema(K)
        if Kmin >= 1 || Kmax <= 1
            K_calculated = false
        end
    end

    if !K_calculated
        K .= suggest_K(model,p,T,z)
    end
end

function tp_flash_fast_K0!(K,model,p,T,z)
    return false
end

function pt_flash_x0(model,p,T,n,method = GeneralizedXYFlash(),non_inx = FillArrays.Fill(false,length(model)),non_iny = FillArrays.Fill(false,length(model));k0 = :wilson)
    ∑n = sum(n)
    z = n/∑n
    if is_vle(method)
        phasex,phasey = :liquid,:vapour
    elseif is_lle(method)
        phasex,phasey = :liquid,:liquid
    else
        phasex,phasey = :unknown,:unknown
    end
    phases = (phasex,phasey)
    non_inw = (non_inx,non_iny)
    nc = length(model)
    _1 = oneunit(Base.promote_eltype(model,p,T,z))
    x,y = fill(_1,nc),fill(_1,nc)
    x .= z
    y .= z
    if !isnothing(method.K0)
        K = _1 * method.K0
        lnK = log.(K)
        volx = zero(_1)
        voly = zero(_1)
    elseif !isnothing(method.x0) && !isnothing(method.y0)
        x = method.x0 ./ sum(method.x0)
        y = method.y0 ./ sum(method.y0)
        lnK = log.(x ./ y)
        volx = zero(_1)
        voly = zero(_1)
        if method.v0 == nothing
            lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,z,nothing,(nothing,nothing),phases,non_inw)
        else
            vl0,vv0 = method.v0
            lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,z,nothing,(vl0,vv0),phases,non_inw)
        end
        K = exp.(lnK)
    elseif is_vle(method) || is_unknown(method)
        # Wilson Correlation for K
        K = tp_flash_K0(model,p,T,z)
        #if we can't predict K, we use lle
        if is_unknown(method)
            Kmin,Kmax = extrema(K)
            if Kmin >= 1 || Kmax <= 1
                K = K0_lle_init(model,p,T,z)
            end
        end
        lnK = log.(K)
        volx = zero(_1)
        voly = zero(_1)
    else
        K = K0_lle_init(model,p,T,z)
        lnK = log.(K)
        volx = zero(_1)
        voly = zero(_1)
    end
    β,status,_ = rachfordrice_β0(K,z,nothing,non_inx,non_iny)
    #if status != RREq, maybe initial K values overshoot the actual phase split.
    if status != RREq
        Kmin,Kmax = K_extrema(K,non_inx,non_iny)
        if !(Kmin >= 1 || Kmax <= 1)
            #valid K, still single phase.
            if status == RRLiquid #bubble point.
                β = eps(typeof(β))
                status = RREq
            elseif status == RRVapour #dew point
                β = one(β) - eps(typeof(β))
                status = RREq
            end
        end
    else
        β = rachfordrice(K, z; β0=β, non_inx=non_inx, non_iny=non_iny)
    end
    if β > 1
        β = one(β)
    elseif β < 0
        β = zero(β)
    end
    y = rr_flash_vapor!(y,K,z,β)
    y ./= sum(y)
    x = rr_flash_liquid!(x,K,z,β)
    x ./= sum(x)
    βv = ∑n*β
    βl = ∑n - βv
    if !isnothing(method.v0) && iszero(volx) && iszero(voly)
        vl0,vv0 = method.v0
        volx,voly = _1*vl0,_1*vv0
    end
    iszero(volx) && (volx = volume(model,p,T,x,phase = phasex))
    iszero(voly) && (voly = volume(model,p,T,y,phase = phasey))
    has_a_res(model) && is_liquid(VT_identify_phase(model,voly,T,y)) && (voly = Rgas(model)*T/p)
    r = FlashResult(p,T,SA[x,y],SA[βl,βv],SA[volx,voly],sort = false)
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