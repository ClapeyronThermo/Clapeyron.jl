"""
    MichelsenTPFlash{T}(;kwargs...)

Method to solve non-reactive multicomponent flash problem by Michelsen's method.

Only two phases are supported. if `K0` is `nothing`, it will be calculated via the Wilson correlation.

### Keyword Arguments:
- `equilibrium`: `:vle` for liquid vapor equilibria, `:lle` for liquid liquid equilibria, `:unknown` if not specified.
- `K0`: initial guess for the K-values.
- `x0`: initial guess for the composition of phase x.
- `y0`: initial guess for the composition of phase y.
- `vol0`: initial guesses for phase x and phase y volumes.
- `K_tol`: tolerance to stop the calculation.
- `ss_iters`: number of Successive Substitution iterations to perform.
- `nacc`: accelerate successive substitution method every nacc steps. Should be a integer bigger than 3. Set to 0 for no acceleration.
- `second_order`: whether to solve the Gibbs energy minimization using the analytical hessian or not.
- `noncondensables`: arrays with names (strings) of components non allowed on the liquid phase. In the case of LLE equilibria, corresponds to the `x` phase.
- `nonvolatiles`: arrays with names (strings) of components non allowed on the vapour phase. In the case of LLE equilibria, corresponds to the `y` phase.
- `flash_result::FlashResult`: can be provided instead of `x0`,`y0` and `vol0` for initial guesses.
"""
struct MichelsenTPFlash{T} <: TPFlashMethod
    equilibrium::Symbol
    K0::Union{Vector{T},Nothing}
    x0::Union{Vector{T},Nothing}
    y0::Union{Vector{T},Nothing}
    v0::Union{Tuple{T,T},Nothing}
    K_tol::Float64
    ss_iters::Int
    nacc::Int
    second_order::Bool
    noncondensables::Union{Nothing,Vector{String}}
    nonvolatiles::Union{Nothing,Vector{String}}
    verbose::Bool
end

michelsen_use_opt_solver(::MichelsenTPFlash) = true
michelsen_itss(method::MichelsenTPFlash) = method.ss_iters

Base.eltype(method::MichelsenTPFlash{T}) where T = T

function index_reduction(m::MichelsenTPFlash,idx::AbstractVector)
    equilibrium,K0,x0,y0,v0,K_tol,ss_iters,nacc,second_order,noncondensables,nonvolatiles,verbose = m.equilibrium,m.K0,m.x0,m.y0,m.v0,m.K_tol,m.ss_iters,m.nacc,m.second_order,m.noncondensables,m.nonvolatiles,m.verbose
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return MichelsenTPFlash(;equilibrium,K0,x0,y0,v0,K_tol,ss_iters,nacc,second_order,noncondensables,nonvolatiles,verbose)
end

numphases(::MichelsenTPFlash) = 2

function MichelsenTPFlash(;equilibrium = :unknown,
                        K0 = nothing,
                        x0 = nothing,
                        y0 = nothing,
                        v0 = nothing,
                        K_tol = sqrt(eps(Float64)),
                        ss_iters = 21,
                        nacc = 5,
                        second_order = false,
                        noncondensables = nothing,
                        nonvolatiles = nothing,
                        flash_result = nothing,
                        verbose = false)
    !(is_vle(equilibrium) | is_lle(equilibrium) | is_unknown(equilibrium))  && throw(error("invalid equilibrium specification for MichelsenTPFlash"))

    if flash_result isa FlashResult
        comps,β,volumes = flash_result.compositions,flash_result.fractions,flash_result.volumes
        np = numphases(flash_result)
        np != 2 && incorrect_np_flash_error(MichelsenTPFlash,flash_result)
        w1,w2 = comps[1],comps[2]
        v = (volumes[1],volumes[2])
        return MichelsenTPFlash(;equilibrium,x0 = w1,y0 = w2,v0 = v,K_tol,ss_iters,nacc,second_order,noncondensables,nonvolatiles,verbose)
    end

    if K0 == x0 == y0 == nothing #nothing specified
        #is_lle(equilibrium)
        T = Nothing
    else
        if !isnothing(K0) & isnothing(x0) & isnothing(y0) #K0 specified
            T = eltype(K0)
        elseif isnothing(K0) & !isnothing(x0) & !isnothing(y0)  #x0, y0 specified
            T = eltype(x0)
        else
            throw(error("invalid specification of initial points"))
        end
    end
    #check for non-volatiles / non-condensables here
    #if is_lle(equilibrium)
    #    if !isnothing(nonvolatiles) && length(nonvolatiles) > 0
    #        throw(error("LLE equilibria does not support setting nonvolatiles"))
    #    end
    #
    #    if !isnothing(noncondensables) && length(noncondensables) > 0
    #        throw(error("LLE equilibria does not support setting noncondensables"))
    #    end
    #end
    nonvolatiles isa String && (nonvolatiles = [nonvolatiles])
    noncondensables isa String && (noncondensables = [noncondensables])
    #check for nacc
    if nacc in (1,2,3) || nacc < 0
        throw(error("incorrect specification for nacc"))
    end

    if T == Nothing && v0 != nothing
        TT = Base.promote_eltype(v0[1],v0[2])
        _v0 = (v0[1],v0[2])
    elseif T != nothing && v0 != nothing
        TT = Base.promote_eltype(one(T),v0[1],v0[2])
        _v0 = (v0[1],v0[2])
    else
        TT = T
        _v0 = v0
    end

    return MichelsenTPFlash{TT}(equilibrium,K0,x0,y0,_v0,K_tol,ss_iters,nacc,second_order,noncondensables,nonvolatiles,verbose)
end

#hook to precalculate things with the activity model.
__tpflash_cache_model(model::EoSModel,p,T,z,equilibrium) = model

__tpflash_gibbs_reduced(model,p,T,x,y,β,eq) = __tpflash_gibbs_reduced(model,p,T,x,y,β,eq,nothing)

function __tpflash_gibbs_reduced(model,p,T,x,y,β,eq,volumes)
    RT = Rgas(model)*T
    if volumes == nothing
        isone(β) && return gibbs_free_energy(model,p,T,y)/RT
        iszero(β) && return gibbs_free_energy(model,p,T,x)/RT
        return (gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/RT
    else
        vx,vy = volumes
        isone(β) && return VT_gibbs_free_energy(model,vy,T,y,p)/RT
        iszero(β) && return VT_gibbs_free_energy(model,vx,T,x,p)/RT
        return (VT_gibbs_free_energy(model,vx,T,x,p)*(1-β)+VT_gibbs_free_energy(model,vy,T,y,p)*β)/RT
    end
end


function tp_flash_impl(model::EoSModel,p,T,z,method::MichelsenTPFlash)

    model_cached = __tpflash_cache_model(model,p,T,z,method.equilibrium)

    x,y,β,v = tp_flash_michelsen(model_cached,p,T,z,method,true)

    volumes = [v[1],v[2]]
    comps = [x,y]
    βi = [1-β ,β]

    if isnan(β)
        g = β
    elseif has_a_res(model_cached)
        g = __tpflash_gibbs_reduced(model_cached,p,T,x,y,β,method.equilibrium,volumes)
    else
        g = __tpflash_gibbs_reduced(model_cached,p,T,x,y,β,method.equilibrium)
    end

    return FlashResult(comps,βi,volumes,FlashData(p,T,g))
end
rrrrr = Ref{Any}()
function tp_flash_michelsen(model::EoSModel, p, T, z, method = MichelsenTPFlash(), reduced = false)

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
    x = similar(z,Base.promote_eltype(model,p,T,z))
    y = similar(z,Base.promote_eltype(model,p,T,z))
    x .= z
    y .= z
    K,lnK = similar(x),similar(x)
    dlnϕ_cache = ∂lnϕ_cache(model, p, T, x, Val{false}())
    _1 = one(eltype(K))
    if !isnothing(K0)
        K .= K0
        lnK .= log.(K)
        verbose && @info "K0 already provided"
    elseif !isnothing(x0) && !isnothing(y0)
        x .= x0 ./ sum(x0)
        y .= y0 ./ sum(y0)
        lnK .= log.(y ./ x)
        lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,z,nothing,(volx,voly),phases,non_inw,dlnϕ_cache)
        K .= exp.(lnK)
        verbose && @info "x0,y0 provided, calculating K0 via Clapeyron.update_K!"
    elseif is_vle(equilibrium) || is_unknown(equilibrium)
        # VLE correlation for K
        verbose && @info "K0 calculated via pure VLE correlation"
        tp_flash_K0!(K,model,p,T,z)

        #if we can't predict K, we use lle
        if is_unknown(equilibrium)
            Kmin,Kmax = K_extrema(K,non_inx,non_iny)
            if Kmin > 1 || Kmax < 1
                verbose && @info "VLE correlation failed, trying LLE initial point."
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
    verbose && @info "K0 = $K"
    Ktype = eltype(K)
    _0 = zero(Ktype)
    _1 = one(Ktype)
    # Initial guess for phase split
    status = rachfordrice_status(K,z,non_inx,non_iny,K_tol = K_tol)
    status0 = status
    #=
    TREND bubble/dew initialization
    Maybe initial K values overshoot the actual phase split.
    if the initial K values generate a single phase result, but we can split the K into two compositions (Kmin < 1 or Kmax > 1)
    then we start at the bubble (or dew conditions)
    =#
    
    if status == RRLiquid
        β = _0
        if maximum(K) >= 1 #liquid phase, but there is posibility to generate a vapour composition
            verbose && @info "suppossing β = 0 (bubble initialization)"
            status = RREq
            β += eps(eltype(β))
        end
    elseif status == RRVapour
        β = _1
        if minimum(K) <= 1 #vapour phase, but there is posibility to generate a liquid composition
            verbose && @info "suppossing β = 1 (dew initialization)"
            status = RREq
            β -= eps(eltype(β))
        end
    elseif status == RREq
        β = rachfordrice(K, z; non_inx=non_inx, non_iny=non_iny,K_tol = K_tol)
    else
        β = _0/_0
    end

    verbose && @info "initial vapour fraction = $β"
    
    verbose && status != RREq && @info "initial point is single-phase (does not satisfy Rachford-Rice constraints). Exiting early"
    # Stage 1: Successive Substitution
    error_lnK = _1
    it = 0
    itacc = 0

    if nacc != 0
        lnK3,lnK4,lnK5,K_dem,lnK_dem,ΔlnK1,ΔlnK2 = similar(lnK),similar(lnK),similar(lnK),similar(lnK),similar(lnK),similar(lnK),similar(lnK)
        x_dem,y_dem = similar(x),similar(y)
    else
        lnK3,lnK4,lnK5,K_dem,lnK_dem,ΔlnK1,ΔlnK2 = lnK,lnK,lnK,lnK,lnK,lnK,lnK
        x_dem,y_dem = x,y
    end

    lnK_old = similar(lnK)
    gibbs = one(_1)
    gibbs_dem = one(_1)
    vcache = Ref((_1, _1))
    verbose && @info "iter         β      error_lnK            K"
    while error_lnK > K_tol && it < itss && status == RREq
        it += 1
        itacc += 1
        lnK_old .= lnK

        x,y = update_rr!(K,β,z,x,y,non_inx,non_iny)

        # Updating K's
        lnK,volx,voly,gibbs = update_K!(lnK,model,p,T,x,y,z,β,(volx,voly),phases,non_inw,dlnϕ_cache)
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
            β_dem = rachfordrice(K_dem, z; β0=β, non_inx=non_inx, non_iny=non_iny)
            x_dem,y_dem = update_rr!(K_dem,β_dem,z,x_dem,y_dem,non_inx,non_iny)
            lnK_dem,volx_dem,voly_dem,gibbs_dem = update_K!(lnK_dem,model,p,T,x_dem,y_dem,z,β_dem,(volx,voly),phases,non_inw,dlnϕ_cache)
            # only accelerate if the Gibbs energy is reduced
            if gibbs_dem < gibbs
                lnK .= lnK_dem
                volx = _1 * volx_dem
                voly = _1 * voly_dem
                vcache[] = (volx,voly)
                β = _1 * β_dem
            end
        end
        K .= exp.(lnK)
        β = rachfordrice(K, z; β0=β, non_inx=non_inx, non_iny=non_iny, K_tol = K_tol)
        verbose && @info "$it    $β  $(round(error_lnK,sigdigits=4)) $K"

        if isnan(β) #try to save K? basically damping
            if rachfordrice_status(K,z,non_inx,non_iny;K_tol = K_tol) != RRTrivial
                K .= 0.5 * K .+ 0.5 * y ./ x
                β = rachfordrice(K, z; non_inx=non_inx, non_iny=non_iny, K_tol = K_tol)
            end
        end

        status = rachfordrice_status(K,z,non_inx,non_iny;K_tol = K_tol)

        Kmin,Kmax = K_extrema(K,non_inx,non_iny)
        if status == RRLiquid && Kmin < 1
            status = RREq
            β = eps(Ktype)
        elseif status == RRVapour && Kmax > 1
            status = RREq
            β = 1 - eps(Ktype)
        end
        # Computing error
        # error_lnK = sum((lnK .- lnK_old).^2)
        error_lnK = dnorm(@view(lnK[in_equilibria]),@view(lnK_old[in_equilibria]),1)
    end
    
    verbose && it > 0 && @info "$it SS iterations done, error(lnK) = $error_lnK"

    if it > 0 && !isnan(β)
        β_margin = min(β, _1 - β)
        verbose && @info "boundary check: β_margin = $(β_margin)"
        # 4 guards to pinpoint the scenario of #465
        if status == RREq && β_margin <= cbrt(K_tol) && Kmin < _1 && Kmax > _1
            status, newβ = rr_margin_check(K,z,non_inx,non_iny;K_tol = K_tol,verbose = verbose)
            status != RREq && (β = newβ)
        end
        # single composition update with the (possibly projected) β
        x, y = update_rr!(K, β, z, x, y, non_inx, non_iny)
    end
    
    # Stage 2: Minimization of Gibbs energy
    if error_lnK > K_tol && it == itss && status == RREq && use_opt_solver
        verbose && @info "$error(lnK) > $K_tol, solving via non-linear system"
        
        nx = similar(K)
        ny = similar(K)
        ny_var0 = y[in_equilibria] * β
        update_nxy!(nx,ny,ny_var0,z,non_inx,non_iny)
        in_eq = (in_equilibria,non_inx,non_iny)
        caches = (nx,ny,vcache,dlnϕ_cache,in_eq,phases)
        flash_obj = michelsen_optimization_obj(model,p,T,z,caches)
        ub = similar(ny_var0)
        ub .= @view z[in_equilibria]
        lb = similar(ny_var0)
        lb .= 0
        opt_options = OptimizationOptions(f_abstol = 1e-12,f_reltol = 1e-12,maxiter = 100)
        if second_order
            sol = Solvers.optimize(flash_obj, ny_var0, Solvers.LineSearch(Solvers.Newton2(ny_var0),Solvers.BoundedLineSearch(lb,ub)),opt_options)
        else
            sol = Solvers.optimize(flash_obj, ny_var0, Solvers.LineSearch(Solvers.BFGS(),Solvers.BoundedLineSearch(lb,ub)),opt_options)
        end

        #= TODO: do something with the values of the optimization procedure
        if abs(sol.info.fx) <= 4*eps(eltype(K))

        elseif sol.info.fx > sol.info.f0 + 4*eps(eltype(K))
        
        end =#
        ny_var = Solvers.x_sol(sol)
        update_nxy!(nx,ny,ny_var,z,non_inx,non_iny)
        x .= nx ./ sum(nx)
        y .= ny ./ sum(ny)
        K .= y ./ x
        β = rachfordrice(K, z; non_inx=non_inx, non_iny=non_iny, K_tol = K_tol)

        β_margin = min(β, _1 - β)
        verbose && @info "boundary check: β_margin = $(β_margin)"
        # 4 guards to pinpoint the scenario of #465
        if β_margin <= cbrt(K_tol) && Kmin < _1 && Kmax > _1
            status, newβ = rr_margin_check(K,z,non_inx,non_iny;K_tol = K_tol,verbose = verbose)
            status != RREq && (β = newβ)
        end
    end
    
    verbose && @info "final K values: $K"
    verbose && @info "final vapour fraction: $β"

    #convergence checks (TODO, seems to fail with activity models)
    status = rachfordrice_status(K,z,non_inx,non_iny,K_tol = K_tol)
    verbose && status != RREq && @info "result is single-phase (does not satisfy Rachford-Rice constraints)."
    vx,vy = vcache[]
    #@show vx,vy
    #maybe azeotrope, do nothing in this case
    if abs(vx - vy) > sqrt(max(abs(vx),abs(vy))) && status != RREq
        verbose && @info "trivial result but different volumes (maybe azeotrope?)"
        status = RREq
    elseif status == RRTrivial && it > 0
        verbose && @info "procedure converged to trivial K-values, checking initial conditions to see if resulting phase is liquid or vapour."
        status0 == RRLiquid && (status = RRLiquid)
        status0 == RRVapour && (status = RRVapour)
    elseif status == RREq && β <= eps(eltype(β))
        status = RRLiquid
    elseif status == RREq && β >=  one(β) - eps(eltype(β))
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

export MichelsenTPFlash
