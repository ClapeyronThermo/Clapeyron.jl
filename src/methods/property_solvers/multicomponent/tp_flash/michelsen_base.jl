function rachfordrice(K, z; β0=nothing, non_inx=FillArrays.Fill(false,length(z)), non_iny=FillArrays.Fill(false,length(z)))
    # Function to solve Rachdord-Rice mass balance
    β,singlephase,limits,_ = rachfordrice_β0(K,z,β0,non_inx,non_iny)
    if length(z) <= 3 && all(Base.Fix2(>,0),z) && all(!,non_inx) && all(!,non_iny) && !singlephase
        return rr_vle_vapor_fraction_exact(K,z)
    end
    #halley refinement
    if !singlephase
        return rr_flash_refine(K,z,β,non_inx,non_iny,limits)
    else
        return β
    end
end

function dgibbs_obj!(model::EoSModel, p, T, z, phasex, phasey,
    nx, ny, vcache, ny_var = nothing, in_equilibria = FillArrays.Fill(true,length(z)), non_inx = in_equilibria, non_iny = in_equilibria;
    F=nothing, G=nothing, H=nothing)

    # Objetive Function to minimize the Gibbs Free Energy
    # It computes the Gibbs free energy, its gradient and its hessian
    iv = 0
    for i in eachindex(z)
        if in_equilibria[i]
            iv += 1
            nyi = ny_var[iv]
            ny[i] = nyi
            nx[i] =z[i] - nyi
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
        # Computing Gibbs Energy Hessian
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
        # Computing Gibbs Energy gradient
        i0 = 0
        for i in eachindex(in_equilibria)
            if in_equilibria[i]
                i0 += 1
                G[i0] = ϕy[i] - ϕx[i]
            end
        end
    end

    if F !== nothing
        # Computing Gibbs Energy
        FO = dot(ny,ϕy) + dot(nx,ϕx)
        return FO
    end
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
            !non_inx[i] || isinf(lnK[i]) && (gibbs += (1-β)*x[i]*(log(x[i]) + lnϕx[i]))
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
            !non_iny[i] || iszero(exp(lnK[i])) && (gibbs += β*y[i]*(log(y[i]) + lnϕy[i]))
        end
    else
        gibbs = gibbs/gibbs
    end
    return lnK,volx,voly,gibbs
end

#updates x,y after a sucessful rachford rice procedure
function update_rr!(K,β,z,x,y,
    non_inx=FillArrays.Fill(false,length(z)),non_iny=non_inx)
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

function tp_flash_K0(model,p,T)
    K = zeros(Base.promote_eltype(model,p,T),length(model))
    return tp_flash_K0!(K,model,p,T)
end

function tp_flash_K0!(K,model,p,T)
    if has_fast_crit_pure(model)
        wilson_k_values!(K,model,p,T)
    else
        pures = split_pure_model(model)
        for i in 1:length(model)
            sat_x = extended_saturation_pressure(pures[i],T)
            K[i] = sat_x[1]/p
        end
    end
    return K
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
    elseif is_vle(method) || is_unknown(method) && k0 == :wilson
        # Wilson Correlation for K
        K = tp_flash_K0(model,p,T)
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
       # volx,voly = NaN*_1,NaN*_1
    elseif is_vle(method) || is_unknown(method)
        K = suggest_K(model,p,T,z)
        lnK = log.(K)
        volx = zero(_1)
        voly = zero(_1)
    else
        K = K0_lle_init(model,p,T,z)
        lnK = log.(K)
        volx = zero(_1)
        voly = zero(_1)
    end
    β,singlephase,_,g01 = rachfordrice_β0(K,z,nothing,non_inx,non_iny)
    g0,g1 = g01
    #if singlephase == true, maybe initial K values overshoot the actual phase split.
    if singlephase
        Kmin,Kmax = extrema(K)
        if !(Kmin >= 1 || Kmax <= 1)
            #valid K, still single phase.
            if g0 <= 0 && g1 < 0 #bubble point.
                β = eps(typeof(β))
                singlephase = false
            elseif g0 > 0 && g1 >= 0 #dew point
                β = one(β) - eps(typeof(β))
                singlephase = false
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
