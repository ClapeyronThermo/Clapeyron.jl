function rachfordrice(K, z; β0=nothing, non_inx=FillArrays.Fill(false,length(z)), non_iny=non_inx)
    # Function to solve Rachdord-Rice mass balance
    β,singlephase,limits = rachfordrice_β0(K,z,β0)
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
function update_K!(lnK,model,p,T,x,y,volx,voly,phasex,phasey,β = nothing,inx = FillArrays.Fill(true,length(x)),iny = inx)
    lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
    lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
    lnK .= lnϕx - lnϕy
    gibbs = zero(eltype(lnK))
    if β !== nothing
        for i in eachindex(y)
            if iny[i]
                gibbs += β*y[i]*(log(y[i]) + lnϕy[i])
            end
            if inx[i]
                gibbs += (1-β)*x[i]*(log(x[i]) + lnϕx[i])
            end
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
        if non_iny[i]
            x[i] = z[i] / (1. - β)
            y[i] = 0.
        end
        # modification for non-in-x components Ki -> ∞
        if non_inx[i]
            x[i] = 0.
            y[i] = z[i] / β
        end
    end
    x ./= sum(x)
    y ./= sum(y)
    return x,y
end

function tp_flash_K0(model,p,T)
    if getparam(model,:Tc) !== nothing
        K = zeros(typeof(p+T+one(eltype(model))),length(model))
        return tp_flash_K0!(K,model,p,T)
    else
        sat_x = initial_points_bd_T.(split_model(model),T)
        psat = first.(sat_x)
        return psat ./ p
    end
end

tp_flash_K0!(K,model,p,T) = wilson_k_values!(K,model,p,T)
