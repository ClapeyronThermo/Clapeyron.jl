function tp_flash_michelsen(model::ElectrolyteModel, p, T, z, method = MichelsenTPFlash(),reduced = false)

    equilibrium = method.equilibrium
    K0 = method.K0
    x0 = method.x0
    y0 = method.y0
    vol0 = method.v0
    K_tol = method.K_tol
    itss = michelsen_itss(method)
    nacc = method.nacc
    second_order = hasfield(typeof(method),:second_order) ? method.second_order : false
    use_opt_solver = michelsen_use_opt_solver(method)
    verbose = method.verbose
    non_inx_list = method.noncondensables
    non_iny_list = method.nonvolatiles

    Z = model.charge
    model_components = component_list(model)
    ions = model_components[Z.!=0]

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
    model_components = component_list(model)
    non_inx = comps_in_equilibria(model_components,non_inx_list)
    non_inx .= (!).(non_inx)
    # constructing non-in-y list
    non_iny = comps_in_equilibria(model_components,non_iny_list)
    non_iny .= (!).(non_iny)

    non_inw = (non_inx,non_iny)
    phases = (phasex,phasey)

    # components that are allowed to be in two phases
    in_equilibria = @. !non_inx & !non_iny

    # Computing the initial guess for the K vector
    TT = Base.promote_eltype(model,p,T,z)
    x = similar(z,TT)
    y = similar(z,TT)
    x .= z
    y .= z
    K,lnK = similar(x),similar(x)
    K̄,lnK̄ = similar(x),similar(x)
    dlnϕ_cache = ∂lnϕ_cache(model, p, T, x, Val{false}())
    _0,_1 = one(TT),one(TT)
    if !isnothing(K0)
        K .= K0
        lnK .= log.(K)
        verbose && @info "K0 already provided"
    elseif !isnothing(x0) && !isnothing(y0)
        x = x0 ./ sum(x0)
        y = y0 ./ sum(y0)
        lnK .= log.(y ./ x)
        lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,z,nothing,(volx,voly),phases,non_inw,dlnϕ_cache)
        K .= exp.(lnK)
        verbose && @info "x0,y0 provided, calculating K0 via Clapeyron.update_K!"
    elseif is_vle(equilibrium) || is_unknown(equilibrium)
        # VLE Correlation for K
        verbose && @info "K0 calculated via pure VLE correlation"
        tp_flash_K0!(K,model,p,T,z)
        #if we can't predict K, we use lle
        if is_unknown(equilibrium)
            Kmin,Kmax = extrema(K)
            if Kmin > 1 || Kmax < 1
                verbose && @info "VLE correlation falied, trying LLE initial point."
                K = K0_lle_init(model,p,T,z)
            end
        end
        lnK .= log.(K)
       # volx,voly = NaN*_1,NaN*_1
    else
        verbose && @info "K0 calculated via LLE initial point (tpd)"
        K .= K0_lle_init(model,p,T,z)
        lnK .= log.(K)
    end
    _1 = one(eltype(K))
    # Initial guess for phase split
    ψ = -sum(Z.*lnK)/sum(abs.(Z))

    #=
    notation:
    K = ϕx/ϕy
    K̄ = K.*exp.(Z.*ψ) = y/x
    =#

    K̄ .= K.*exp.(Z.*ψ)
    logK̄ .= log.(K̄)
    β,status,_ = rachfordrice_β0(K̄,z,nothing,non_inx,non_iny)
    status0 = status
    #=TODO:
    there is a method used in TREND that tries to obtain adequate values of K
    in the case of incorrect initialization.
    =#
    # Stage 1: Successive Substitution
    verbose && @info "initial vapour fraction = $β"
    verbose && @info "ψ(K0) = $ψ"
    verbose && status != RREq && @info "initial point is single-phase (does not satisfy Rachford-Rice constraints). Exiting early"

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

    lnK̄_old = similar(lnK)
    β_old = typemax(TT)
    gibbs = one(_1)
    gibbs_dem = one(_1)
    vcache = Ref((_1, _1))
    verbose && @info "iter  status        β      error_lnK            K"
    while (error_lnK > K_tol || abs(β_old-β) > 1e-9) && it < itss && status in (RREq,RRLiquid,RRVapour)
        it += 1
        itacc += 1
        lnK̄_old .= lnK + Z.*ψ
        β_old = β

        x,y = update_rr!(K̄,β,z,x,y,non_inx,non_iny)

        # Updating K's
        lnK,volx,voly,gibbs = update_K!(lnK,model,p,T,x,y,z,β,(volx,voly),phases,non_inw,dlnϕ_cache)
        vcache[] = (volx,voly)
        gibbs +=  β*ψ*dot(x,Z)

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
            β_dem,ψ_dem = rachfordrice(K_dem, z, Z; β0=β, ψ0=ψ, non_inx=non_inx, non_iny=non_iny, verbose = verbose)
            K̄_dem .= K_dem .* exp.(Z .* ψ_dem)
            x_dem,y_dem = update_rr!(K̄_dem,β_dem,z,x_dem,y_dem,non_inx,non_iny)
            lnK_dem,volx_dem,voly_dem,gibbs_dem = update_K!(lnK_dem,model,p,T,x_dem,y_dem,z,β_dem,(volx,voly),phases,non_inw,dlnϕ_cache)
            #add effect of electroneutrality condition on Gibbs energy
            gibbs_dem += β_dem*ψ_dem*dot(x_dem,Z)
            # only accelerate if the Gibbs energy is reduced
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
        β,ψ = rachfordrice(K, z, Z; β0=β, ψ0=ψ, non_inx=non_inx, non_iny=non_iny,verbose = verbose)
        

        lnK̄ .= lnK .+ Z.*ψ
        K̄ = exp.(lnK̄)
        status = rachfordrice_status(K̄,z,non_inx,non_iny;K_tol)
    
        verbose && @info "$it    $status   $β  $(round(error_lnK,sigdigits=4)) $K̄"
    
        # Computing error
        # error_lnK = sum((lnK .- lnK_old).^2)
        error_lnK = dnorm(@view(lnK̄[in_equilibria]),@view(lnK̄_old[in_equilibria]),1)
    end

    verbose && it > 0 && @info "$it SS iterations done, error(lnK) = $error_lnK"

    if it > 0 && !isnan(β)
        # single composition update with the (possibly projected) β
        x, y = update_rr!(K̄, β, z, x, y, non_inx, non_iny)
    end

    # Stage 2: Minimization of Gibbs energy
    if error_lnK > K_tol && it == itss && status == RREq && use_opt_solver
        verbose && @info "$it error(lnK) > $K_tol, solving via non-linear system"
        nx = similar(K)
        ny = similar(K)
        ny_var_and_ψ0 = similar(K,count(in_equilibria)+1)
        ny_var_and_ψ0[1:end-1] .= @view(y[in_equilibria]) .* β
        ny_var_and_ψ0[end] = ψ
        update_nxy!(nx,ny,@view(ny_var_and_ψ0[1:end-1]),z,non_inx,non_iny)
        in_eq = (in_equilibria,non_inx,non_iny)
        caches = (nx,ny,vcache,dlnϕ_cache,in_eq,phases)
        flash_obj = michelsen_optimization_obj(model,p,T,z,caches)
        ub = similar(ny_var_and_ψ0)
        ub[1:end-1] .= @view z[in_equilibria]
        lb = similar(ny_var_and_ψ0)
        lb .= 0
        lb[end] = -Inf
        ub[end] = Inf
        opt_options = OptimizationOptions(f_abstol = 1e-12,f_reltol = 1e-8,maxiter = 100)
        if second_order
            sol = Solvers.optimize(flash_obj, ny_var_and_ψ0, Solvers.LineSearch(Solvers.Newton2(ny_var0),Solvers.BoundedLineSearch(lb,ub)),opt_options)
        else
            sol = Solvers.optimize(flash_obj, ny_var_and_ψ0, Solvers.LineSearch(Solvers.BFGS(),Solvers.BoundedLineSearch(lb,ub,)),opt_options)
        end
        ny_var_and_ψ = Solvers.x_sol(sol)
        ny_var = @view ny_var_and_ψ[1:end-1]
        ψ = ny_var_and_ψ[end]
        update_nxy!(nx,ny,ny_var,z,non_inx,non_iny)
        x .= nx ./ sum(nx)
        y .= ny ./ sum(ny)
        K̄ .= y ./ x
        β = rachfordrice(K̄, z; non_inx, non_iny, K_tol, verbose)
    end

    verbose && @info "final K values: $K"

    verbose && @info "final vapour fraction: $β"
    #convergence checks (TODO, seems to fail with activity models)
    status = rachfordrice_status(K̄,z,non_inx,non_iny;K_tol = K_tol)
    verbose && status != RREq && @info "result is single-phase (does not satisfy Rachford-Rice constraints)."

    vx,vy = vcache[]
    #@show vx,vy
    #maybe azeotrope, do nothing in this case
    if abs(vx - vy) > sqrt(max(abs(vx),abs(vy))) && status != RREq
        verbose && @info "trivial result but different volumes (maybe azeotrope?)"
        status = RREq
    elseif status == RRTrivial
        verbose && @info "procedure converged to trivial K-values, checking initial conditions to see if resulting phase is liquid or vapour."
        status0 == RRLiquid && (status = RRLiquid)
        status0 == RRVapour && (status = RRVapour)
    elseif status == RREq && β <= eps(eltype(β))
        status = RRLiquid
    elseif status == RREq && β >= one(β) - eps(eltype(β))
        status = RRVapour
    elseif !material_balance_rr_converged((x,y),z,β) #material balance failed
        verbose && @info "material balance failed."
        status = RRFailure
    end

    verbose && status == RRLiquid && @info "procedure converged to a single liquid phase."
    verbose && status == RRVapour && @info "procedure converged to a single vapour phase."

    if status != RREq
        _0 = zero(eltype(x))
        _1 = one(eltype(x))
        x .= z
        y .= z
        if status == RRLiquid
            β = _0
            vz = volume(model,p,T,z,phase = :l)
        elseif status == RRVapour
            β = _1
            vz = volume(model,p,T,z,phase = :v)
        else
            β = _0/_0
            vz = _0/_0
        end
        vx = vz
        vy = vz
    end

    if !reduced
        x = index_expansion(x,z_nonzero)
        y = index_expansion(y,z_nonzero)
    end
    return x, y, β, (vx,vy)
end

function rachfordrice(K, z, Z; β0=nothing, ψ0=nothing, non_inx=FillArrays.Fill(false,length(z)), non_iny=FillArrays.Fill(false,length(z)),verbose = false)
    # Function to solve Rachdord-Rice mass balance
    β,status,limits = rachfordrice_β0(K.*exp.(Z.*ψ0),z,β0,non_inx,non_iny)
    if status == RREq
        function rachford_rice_donnan(x,K,z,Z)
            β = x[1]
            ψ = x[2]
            F1 = zero(Base.promote_eltype(K,z))
            F2 = zero(Base.promote_eltype(K,z))
            for i in 1:length(Z)
                Zi,zi = Z[i],z[i]
                K̄i =  K[i]*exp(Zi*ψ)
                F1 += zi*(1 - K̄i)/(1 + β*(K̄i - 1)) #rachford rice
                if Zi != 0
                    F2 += zi*Zi/(1 + β*(K̄i - 1)) #electroneutrality of phase x
                end
            end
            return SVector((F1,F2))
        end
        x0 = SVector(Base.promote(β0,ψ0))
        ff(F,x) = rachford_rice_donnan(x,K,z,Z)
        results = Solvers.nlsolve(ff,x0)
        sol = Clapeyron.Solvers.x_sol(results)
        β = sol[1]
        ψ = sol[2]
        return SVector(Base.promote(β,ψ))
    else
        return SVector(Base.promote(β0,ψ0))
    end
end

function michelsen_optimization_of!(g,H,model::ElectrolyteModel,p,T,z,caches,ny_var_and_ψ,gz)
    ny_var = @view ny_var_and_ψ[1:end-1]
    ψ = ny_var[end]

    second_order = !isnothing(H)
    nx,ny,vcache,lnϕ_cache,in_eq,phases = caches
    in_equilibria,non_inx,non_iny = in_eq
    phasex,phasey = phases
    volx,voly = vcache[]
    iv = 0
    Z = model.charge
    update_nxy!(nx,ny,ny_var,z,non_inx,non_iny)
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

        !isnothing(g) && (g[1:end-1] .= @view ∂x[in_equilibria])
        !isnothing(g) && (g[1:end-1] .+= @view(Z[in_equilibria]) .* ψ)

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

        !isnothing(g) && (g[1:end-1] .-= @view ∂y[in_equilibria])
        !isnothing(g) && (g[1:end-1] .*= -1)

        f += dot(∂y,ny) + dot(nx,Z)
    else
        ∂x,volx = lnϕ(model, p, T, x,lnϕ_cache; phase=phasex, vol0=volx)
        for i in 1:size(∂x,1)
            ∂x[i] += log(x[i])
            non_inx[i] && (∂x[i] = 0)
        end
        !isnothing(g) && (g[1:end-1] .= @view ∂x[in_equilibria])
        !isnothing(g) && (g[1:end-1] .+= @view(Z[in_equilibria]) .* ψ)
        f += dot(∂x,nx)
        ∂y,voly = lnϕ(model, p, T, y,lnϕ_cache; phase=phasey, vol0=voly)
        for i in 1:size(∂y,1)
            ∂y[i] += log(y[i])
            non_iny[i] && (∂y[i] = 0)
        end
        !isnothing(g) && (g[1:end-1] .-= @view ∂y[in_equilibria])
        !isnothing(g) && (g[1:end-1] .*= -1)

        f += dot(∂y,ny) + dot(nx,Z)
    end
    vcache[] = (volx,voly)
    return f
end

function michelsen_gibbs_feed(model::ElectrolyteModel,p,T,z,caches)
    nx,ny,vcache,lnϕ_cache,in_eq,phases = caches
    in_equilibria,non_inx,non_iny = in_eq
    ∂z,volz = lnϕ(model, p, T, z, lnϕ_cache,phase = :l)
    ∂z  .+= log.(z)
    gz = dot(@view(z[in_equilibria]),@view(∂z[in_equilibria]))
    return gz
end