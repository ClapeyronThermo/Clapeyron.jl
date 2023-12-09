function tp_flash_michelsen(model::ElectrolyteModel, p, T, z; equilibrium=:vle, K0=nothing,
                                     x0=nothing, y0=nothing, vol0=(nothing, nothing),
                                     K_tol=1e-12, itss=101, nacc=5, second_order=false, use_opt_solver = true,
                                     non_inx_list=nothing, non_iny_list=nothing, reduced=false)


    Z = model.charge
    ions = model.components[Z.!=0]
    
    if !reduced
        model_full,z_full = model,z
        model,z_nonzero = index_reduction(model_full,z_full)
        z = z_full[z_nonzero]
    end

    if is_vle(equilibrium)
        phasex = :liquid
        phasey = :vapor
        append!(non_iny_list,ions)
    elseif is_lle(equilibrium)
        phasex = :liquid
        phasey = :liquid
    end

    # Setting the initial guesses for volumes
    vol0 === nothing && (vol0 = (nothing,nothing))
    volx, voly = vol0

    nc = length(model)
    # constructing non-in-x list
    if !isnothing(non_inx_list)
        non_inx_names_list = [x for x in non_inx_list if x in model.components]
    else
        non_inx_names_list = String[]
    end

    if !isnothing(non_iny_list)
        non_iny_names_list = [x for x in non_iny_list if x in model.components]
    else
        non_iny_names_list = String[]
    end

    # constructing non-in-x list
    non_inx = Bool.(zeros(nc))
    # constructing non-in-y list
    non_iny = Bool.(zeros(nc))

    for i in 1:nc
        component = model.components[i]
        if component in non_inx_names_list
            non_inx[i] = true
        end

        if component in non_iny_names_list
            non_iny[i] = true
        end
    end

    inx = .!non_inx
    iny = .!non_iny

    active_inx = !all(inx)
    active_iny = !all(iny)
    # components that are allowed to be in two phases
    in_equilibria = inx .& iny
    # Computing the initial guess for the K vector
    x = similar(z)
    y = similar(z)
    if !isnothing(K0)
        K = 1. * K0
        lnK = log.(K)
    elseif !isnothing(x0) && !isnothing(y0)
        x = x0 ./ sum(x0)
        y = y0 ./ sum(y0)
        lnK = log.(x ./ y)
        lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,volx,voly,phasex,phasey,nothing,inx,iny)
        K = exp.(lnK)
    elseif is_vle(equilibrium)
        # Wilson Correlation for K
        K = tp_flash_K0(model,p,T)
        lnK = log.(K)
       # volx,voly = NaN*_1,NaN*_1
    else
        K = K0_lle_init(model,p,T,z)
        lnK = log.(K)
    end
    _1 = one(p+T+first(z))
    # Initial guess for phase split
    ψ = 0.
    β,singlephase,_ = rachfordrice_β0(K.*exp.(Z.*ψ),z)
    #=TODO:
    there is a method used in TREND that tries to obtain adequate values of K
    in the case of incorrect initialization.
    =#
    # Stage 1: Successive Substitution
    error_lnK = _1
    it = 0
    x_dem = similar(z)
    y_dem = similar(z)
    itacc = 0
    lnK3 = similar(lnK)
    lnK4 = similar(lnK)
    lnK5 = similar(lnK)
    K_dem = similar(lnK)
    lnK_dem = similar(lnK)
    ΔlnK1 = similar(lnK)
    ΔlnK2 = similar(lnK)
    gibbs = one(_1)
    gibbs_dem = one(_1)
    vcache = Ref((_1, _1))
    while error_lnK > K_tol && it < itss && !singlephase
        it += 1
        itacc += 1
        lnK_old = lnK .* _1
        β,ψ = rachfordrice(K, z, Z; β0=β, ψ0=ψ, non_inx=non_inx, non_iny=non_iny)
        singlephase = !(0 < β < 1) #rachford rice returns 0 or 1 if it is single phase.
        x,y = update_rr!(K.*exp.(Z.*ψ),β,z,x,y,non_inx,non_iny)
        # Updating K's
        lnK,volx,voly,gibbs = update_K!(lnK,model,p,T,x,y,volx,voly,phasex,phasey,β,inx,iny)
        vcache[] = (volx,voly)
        # acceleration step
        if itacc == (nacc - 2)
            lnK3 = 1. * lnK
        elseif itacc == (nacc - 1)
            lnK4 = 1. * lnK
        elseif itacc == nacc
            itacc = 0
            lnK5 = 1. * lnK
            # acceleration using DEM (1 eigenvalues)
            lnK_dem = dem!(lnK_dem, lnK5, lnK4, lnK3,(ΔlnK1,ΔlnK2))
            K_dem .= exp.(lnK_dem)
            β_dem,ψ_dem = rachfordrice(K_dem, z, Z; β0=β, ψ0=ψ, non_inx=non_inx, non_iny=non_iny)
            x_dem,y_dem = update_rr!(K_dem.*exp.(Z.*ψ_dem),β_dem,z,x_dem,y_dem,non_inx,non_iny)
            lnK_dem,volx_dem,voly_dem,gibbs_dem = update_K!(lnK_dem,model,p,T,x_dem,y_dem,volx,voly,phasex,phasey,β,inx,iny)
            # only accelerate if the gibbs free energy is reduced
            if gibbs_dem < gibbs
                lnK .= _1 * lnK_dem
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
        error_lnK = dnorm(lnK,lnK_old,1)
    end
    K .= x ./ y
    #convergence checks (TODO, seems to fail with activity models)
    _,singlephase,_ = rachfordrice_β0(K.*exp.(Z.*ψ),z)
    vx,vy = vcache[]
    #@show vx,vy
    #maybe azeotrope, do nothing in this case
    if abs(vx - vy) > sqrt(max(abs(vx),abs(vy))) && singlephase
        singlephase = false
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

    if vx < vy #sort by increasing volume
        return x, y, β
    else
        return y, x, 1 - β
    end
end

function rachfordrice(K, z, Z; β0=nothing, ψ0=nothing, non_inx=FillArrays.Fill(false,length(z)), non_iny=non_inx)
    # Function to solve Rachdord-Rice mass balance
    β,singlephase,limits = rachfordrice_β0(K.*exp.(Z.*ψ0),z,β0)
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