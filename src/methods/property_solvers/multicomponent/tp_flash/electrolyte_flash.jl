function tp_flash_michelsen(model::ElectrolyteModel, p, T, z; equilibrium=:vle, K0=nothing,
                                     x0=nothing, y0=nothing, vol0=(nothing, nothing),
                                     K_tol=1e-12, itss=21, nacc=5, second_order=false, use_opt_solver = true,
                                     non_inx_list=nothing, non_iny_list=nothing, reduced=false)


    Z = model.charge
    ions = model.components[Z.!=0]
    
    if !reduced
        model_full,z_full = model,z
        model,z_nonzero = index_reduction(model_full,z_full)
        z = z_full[z_nonzero]
    end

    if is_vle(equilibrium) || is_unknown(equilibrium)
        phasex,phasey = :liquid,:vapour
    elseif is_lle(equilibrium)
        phasex,phasey = :liquid,:liquid
    end

    # Setting the initial guesses for volumes
    vol0 === nothing && (vol0 = (nothing,nothing))
    volx, voly = vol0

    nc = length(model)
    # constructing non-in-x list
    non_inx = fill(false,nc)
    # constructing non-in-y list
    non_iny = fill(false,nc)

    for i in 1:nc
        component = model.components[i]
        non_inx[i] = !isnothing(non_inx_list) && (component in non_inx_list) && true
        non_iny[i] = !isnothing(non_iny_list) && (component in non_iny_list) && true
    end

    non_inw = (non_inx,non_iny)
    phases = (phasex,phasey)

    # components that are allowed to be in two phases
    in_equilibria = @. !non_inx & !non_iny
    # Computing the initial guess for the K vector
    x = similar(z)
    y = similar(z)
    dlnϕ_cache = ∂lnϕ_cache(model, p, T, x, Val{false}())
    if !isnothing(K0)
        K = 1. * K0
        lnK = log.(K)
    elseif !isnothing(x0) && !isnothing(y0)
        x = x0 ./ sum(x0)
        y = y0 ./ sum(y0)
        lnK = log.(y ./ x)
        lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,nothing,(volx,voly),phases,non_inw)
        K = exp.(lnK)
    elseif is_vle(equilibrium) || is_unknown(equilibrium)
        # Wilson Correlation for K
        K = tp_flash_K0(model,p,T)
        #if we can't predict K, we use lle
        if is_unknown(equilibrium)
            Kmin,Kmax = extrema(K)
            
            if Kmin > 1 || Kmax < 1 
                K = K0_lle_init(model,p,T,z)
            end
        end
        lnK = log.(K)
       # volx,voly = NaN*_1,NaN*_1
    else
        K = K0_lle_init(model,p,T,z)
        lnK = log.(K)
    end
    _1 = one(p+T+first(z))
    # Initial guess for phase split
    ψ = 0.
    K̄ = K.*exp.(Z.*ψ)
    β,singlephase,_,_ = rachfordrice_β0(K̄,z,nothing,non_inx,non_iny)
    
    #=TODO:
    there is a method used in TREND that tries to obtain adequate values of K
    in the case of incorrect initialization.
    =#
    # Stage 1: Successive Substitution
    error_lnK = _1
    it = 0
    itacc = 0
    
    if nacc != 0
        lnK3,lnK4,lnK5,K_dem,lnK_dem,ΔlnK1,ΔlnK2,K̄_dem = similar(lnK),similar(lnK),similar(lnK),similar(lnK),similar(lnK),similar(lnK),similar(lnK),similar(lnK)
        x_dem,y_dem = similar(x),similar(y)
    else
        lnK3,lnK4,lnK5,K_dem,lnK_dem,ΔlnK1,ΔlnK2,K̄_dem = lnK,lnK,lnK,lnK,lnK,lnK,lnK,lnK
        x_dem,y_dem = x,y
    end

    lnK_old = similar(lnK)
    gibbs = one(_1)
    gibbs_dem = one(_1)
    vcache = Ref((_1, _1))
    while error_lnK > K_tol && it < itss && !singlephase
        it += 1
        itacc += 1
        lnK_old .= lnK .* _1
        β,ψ = rachfordrice(K, z, Z; β0=β, ψ0=ψ, non_inx=non_inx, non_iny=non_iny)
        singlephase = !(0 < β < 1) #rachford rice returns 0 or 1 if it is single phase.
        K̄ .= K.*exp.(Z.*ψ)
        x,y = update_rr!(K̄,β,z,x,y,non_inx,non_iny)
        # Updating K's
        lnK,volx,voly,gibbs = update_K!(lnK,model,p,T,x,y,β,(volx,voly),phases,non_inw,dlnϕ_cache)
        vcache[] = (volx,voly)
        # acceleration step
        if itacc == (nacc - 2)
            lnK3 .= lnK
        elseif itacc == (nacc - 1)
            lnK4 .= lnK
        elseif itacc == nacc
            itacc = 0
            lnK5 .= lnK
            # acceleration using DEM (1 eigenvalues)
            lnK_dem = dem!(lnK_dem, lnK5, lnK4, lnK3,(ΔlnK1,ΔlnK2))
            K_dem .= exp.(lnK_dem)
            β_dem,ψ_dem = rachfordrice(K_dem, z, Z; β0=β, ψ0=ψ, non_inx=non_inx, non_iny=non_iny)
            K̄ .= K_dem.*exp.(Z.*ψ_dem)
            K̄_dem .= K̄
            x_dem,y_dem = update_rr!(K̄_dem,β_dem,z,x_dem,y_dem,non_inx,non_iny)
            lnK_dem,volx_dem,voly_dem,gibbs_dem = update_K!(lnK_dem,model,p,T,x_dem,y_dem,β,(volx,voly),phases,non_inw,dlnϕ_cache)
            # only accelerate if the gibbs free energy is reduced
            if gibbs_dem < gibbs
                lnK .= lnK_dem
                volx = _1 * volx_dem
                voly = _1 * voly_dem
                vcache[] = (volx,voly)
                β = _1 * β_dem
                ψ = _1 * ψ_dem
            end
        end
        K .= exp.(lnK)

        # Computing error
        # error_lnK = sum((lnK .- lnK_old).^2)
        dnorm(@view(lnK[in_equilibria]),@view(lnK_old[in_equilibria]),1)
    end

    if error_lnK > K_tol && it == itss && !singlephase && use_opt_solver
        nx = zeros(nc)
        ny = zeros(nc)

        if any(non_inx)
            ny[non_inx] = z[non_inx]
            nx[non_inx] .= 0.
        end
        if any(non_iny)
            ny[non_iny] .= 0.
            nx[non_iny] = z[non_iny]
        end

        ny_var0 = y[in_equilibria] * β
        fgibbs!(F, G, H, ny_var) = dgibbs_obj!(model, p, T, z, phasex, phasey,
                                                        nx, ny, vcache, ny_var[1:end-1], ny_var[end] , in_equilibria, non_inx, non_iny;
                                                        F=F, G=G, H=H)
        append!(ny_var0,ψ)

        fgibbs!(F, G, ny_var) = fgibbs!(F, G, nothing, ny_var)

        if second_order
            sol = Solvers.optimize(Solvers.only_fgh!(fgibbs!), ny_var0, Solvers.LineSearch(Solvers.Newton()))
        else
            sol = Solvers.optimize(Solvers.only_fg!(fgibbs!), ny_var0, Solvers.LineSearch(Solvers.BFGS()))
        end
        ny_var = Solvers.x_sol(sol)[1:end-1]
        ψ = Solvers.x_sol(sol)[end]
        ny[in_equilibria] = ny_var
        nx[in_equilibria] .= @view(z[in_equilibria]) .- @view(ny[in_equilibria])
        nxsum = sum(nx)
        nysum = sum(ny)
        x .= nx ./ nxsum
        y .= ny ./ nysum
        β = sum(ny)

    end
    K .= y ./ x
    K̄ .= K.*exp.(Z.*ψ)
    #convergence checks (TODO, seems to fail with activity models)
    _,singlephase,_,_ = rachfordrice_β0(K̄,z,β,non_inx,non_iny)
    vx,vy = vcache[]
    #@show vx,vy
    #maybe azeotrope, do nothing in this case
    if abs(vx - vy) > sqrt(max(abs(vx),abs(vy))) && singlephase
        singlephase = false
    elseif any(isnan,view(K,in_equilibria)) || isnan(ψ)
        singlephase = true
    end
    if singlephase
        β = zero(β)/zero(β)
        x .= z
        y .= z
    end

    if !reduced
        x = index_expansion(x,z_nonzero)
        y = index_expansion(y,z_nonzero)
    end
    return x, y, β, (vx,vy)
end

function rachfordrice(K, z, Z; β0=nothing, ψ0=nothing, non_inx=FillArrays.Fill(false,length(z)), non_iny=FillArrays.Fill(false,length(z)))
    # Function to solve Rachdord-Rice mass balance
    β,singlephase,limits,_ = rachfordrice_β0(K.*exp.(Z.*ψ0),z,β0,non_inx,non_iny)
    if !singlephase
        function rachford_rice_donnan(F,x,K,z,Z)
            β = x[1]
            ψ = x[2]
            F[1] = sum((z.*(1 .-K.*exp.(Z*ψ)))./(1 .+β*(K.*exp.(Z*ψ).-1)))
            F[2] = sum((z.*Z./(1 .+β*(K.*exp.(Z*ψ).-1))))
        end
        f!(F,x) = rachford_rice_donnan(F,x,K,z,Z)
        results = Solvers.nlsolve(f!,[β0,ψ0],TrustRegion(Newton(), Dogleg()))
        sol = Clapeyron.Solvers.x_sol(results)
        β = sol[1]
        ψ = sol[2]
        return β, ψ
    else
        return β, ψ0
    end
end

function dgibbs_obj!(model::ElectrolyteModel, p, T, z, phasex, phasey,
    nx, ny, vcache, ny_var = nothing, ψ = nothing, in_equilibria = FillArrays.Fill(true,length(z)), non_inx = in_equilibria, non_iny = in_equilibria;
    F=nothing, G=nothing, H=nothing)

    Z = model.charge
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
                G[i0] = ϕy[i] - ϕx[i] - ψ*Z[i]
            end
        end
        G[i0+1] = dot(ny,Z)
    end

    if F !== nothing
        # Computing Gibbs Energy
        FO = dot(ny,ϕy) + dot(nx,ϕx) + ψ*dot(nx,Z)
        return FO
    end
end
