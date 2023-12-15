"""
    MichelsenTPFlash{T}(;kwargs...)

Method to solve non-reactive multicomponent flash problem by Michelsen's method.

Only two phases are supported. if `K0` is `nothing`, it will be calculated via the Wilson correlation.

### Keyword Arguments:
- equilibrium = equilibrium type ":vle" for liquid vapor equilibria, ":lle" for liquid liquid equilibria
- `K0` (optional), initial guess for the constants K
- `x0` (optional), initial guess for the composition of phase x
- `y0` = optional, initial guess for the composition of phase y
- `vol0` = optional, initial guesses for phase x and phase y volumes
- `K_tol` = tolerance to stop the calculation
- `ss_iters` = number of Successive Substitution iterations to perform
- `nacc` =  accelerate successive substitution method every nacc steps. Should be a integer bigger than 3. Set to 0 for no acceleration.
- `second_order` = wheter to solve the gibbs energy minimization using the analytical hessian or not
- `noncondensables` = arrays with names (strings) of components non allowed on the liquid phase. In the case of LLE equilibria, corresponds to the `x` phase
- `nonvolatiles` = arrays with names (strings) of components non allowed on the vapour phase. In the case of LLE equilibria, corresponds to the `y` phase

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
end

Base.eltype(method::MichelsenTPFlash{T}) where T = T

function index_reduction(m::MichelsenTPFlash,idx::AbstractVector)
    equilibrium,K0,x0,y0,v0,K_tol,ss_iters,nacc,second_order,noncondensables,nonvolatiles = m.equilibrium,m.K0,m.x0,m.y0,m.v0,m.K_tol,m.ss_iters,m.nacc,m.second_order,m.noncondensables,m.nonvolatiles
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return MichelsenTPFlash(;equilibrium,K0,x0,y0,v0,K_tol,ss_iters,nacc,second_order,noncondensables,nonvolatiles)
end

numphases(::MichelsenTPFlash) = 2

function MichelsenTPFlash(;equilibrium = :vle,
                        K0 = nothing,
                        x0 = nothing,
                        y0 = nothing,
                        v0 = nothing,
                        K_tol = sqrt(eps(Float64)),
                        ss_iters = 21,
                        nacc = 5,
                        second_order = false,
                        noncondensables = nothing,
                        nonvolatiles = nothing)
    !(is_vle(equilibrium) | is_lle(equilibrium)) && throw(error("invalid equilibrium specification for MichelsenTPFlash"))
    if K0 == x0 == y0 === v0 == nothing #nothing specified
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

    #check for nacc
    if nacc in (1,2,3) || nacc < 0
        throw(error("incorrect specification for nacc"))
    end

    return MichelsenTPFlash{T}(equilibrium,K0,x0,y0,v0,K_tol,ss_iters,nacc,second_order,noncondensables,nonvolatiles)
end

is_vle(method::MichelsenTPFlash) = is_vle(method.equilibrium)
is_lle(method::MichelsenTPFlash) = is_lle(method.equilibrium)

#hook to precalculate things with the activity model.
__tpflash_cache_model(model::EoSModel,p,T,z,equilibrium) = model

function __tpflash_gibbs_reduced(model,p,T,x,y,β,eq)
    (gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/Rgas(model)/T
end

function tp_flash_impl(model::EoSModel,p,T,z,method::MichelsenTPFlash)

    model_cached = __tpflash_cache_model(model,p,T,z,method.equilibrium)

    x,y,β =  tp_flash_michelsen(model_cached,p,T,z;equilibrium = method.equilibrium, K0 = method.K0,
            x0 = method.x0, y0 = method.y0, vol0 = method.v0,
            K_tol = method.K_tol,itss = method.ss_iters, nacc=method.nacc,
            second_order = method.second_order,
            non_inx_list=method.noncondensables, non_iny_list=method.nonvolatiles,
            reduced = true)

    G = __tpflash_gibbs_reduced(model_cached,p,T,x,y,β,method.equilibrium)
    X = hcat(x,y)'
    nvals = X.*[1-β
                β] .* sum(z)
    return (X, nvals, G)
end


function tp_flash_michelsen(model::EoSModel, p, T, z; equilibrium=:vle, K0=nothing,
                                     x0=nothing, y0=nothing, vol0=(nothing, nothing),
                                     K_tol=1e-8, itss=21, nacc=5, second_order=false, use_opt_solver = true,
                                     non_inx_list=nothing, non_iny_list=nothing, reduced=false)


    if !reduced
        model_full,z_full = model,z
        model,z_nonzero = index_reduction(model_full,z_full)
        z = z_full[z_nonzero]
    end

    if is_vle(equilibrium)
        phasex = :liquid
        phasey = :vapor
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
    β,singlephase,_ = rachfordrice_β0(K,z)
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
        β = rachfordrice(K, z; β0=β, non_inx=non_inx, non_iny=non_iny)
        singlephase = !(0 < β < 1) #rachford rice returns 0 or 1 if it is single phase.
        x,y = update_rr!(K,β,z,x,y,non_inx,non_iny)
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
            β_dem = rachfordrice(K_dem, z; β0=β, non_inx=non_inx, non_iny=non_iny)
            x_dem,y_dem = update_rr!(K_dem,β_dem,z,x_dem,y_dem,non_inx,non_iny)
            lnK_dem,volx_dem,voly_dem,gibbs_dem = update_K!(lnK_dem,model,p,T,x_dem,y_dem,volx,voly,phasex,phasey,β,inx,iny)
            # only accelerate if the gibbs free energy is reduced
            if gibbs_dem < gibbs
                lnK .= _1 * lnK_dem
                volx = _1 * volx_dem
                voly = _1 * voly_dem
                vcache[] = (volx,voly)
                β = _1 * β_dem
            end
        end
        K .= exp.(lnK)

        # Computing error
        # error_lnK = sum((lnK .- lnK_old).^2)
        error_lnK = dnorm(lnK,lnK_old,1)
    end
    # Stage 2: Minimization of Gibbs Free Energy
    if error_lnK > K_tol && it == itss && !singlephase && use_opt_solver
        nx = zeros(nc)
        ny = zeros(nc)

        if active_inx
            ny[non_inx] = z[non_inx]
            nx[non_inx] .= 0.
        end
        if active_iny
            ny[non_iny] .= 0.
            nx[non_iny] = z[non_iny]
        end

        ny_var0 = y[in_equilibria] * β
        fgibbs!(F, G, H, ny_var) = dgibbs_obj!(model, p, T, z, phasex, phasey,
                                                        nx, ny, vcache, ny_var, in_equilibria, non_inx, non_iny;
                                                        F=F, G=G, H=H)

        fgibbs!(F, G, ny_var) = fgibbs!(F, G, nothing, ny_var)

        if second_order
            sol = Solvers.optimize(Solvers.only_fgh!(fgibbs!), ny_var0, Solvers.LineSearch(Solvers.Newton()))
        else
            sol = Solvers.optimize(Solvers.only_fg!(fgibbs!), ny_var0, Solvers.LineSearch(Solvers.BFGS()))
        end
        ny_var = Solvers.x_sol(sol)
        ny[in_equilibria] = ny_var
        nx[in_equilibria] = z[in_equilibria] .- ny[in_equilibria]
        nxsum = sum(nx)
        nysum = sum(ny)
        x = nx ./ nxsum
        y = ny ./ nysum
        β = sum(ny)

    end
    K .= x ./ y
    #convergence checks (TODO, seems to fail with activity models)
    _,singlephase,_ = rachfordrice_β0(K,z)
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

export MichelsenTPFlash
