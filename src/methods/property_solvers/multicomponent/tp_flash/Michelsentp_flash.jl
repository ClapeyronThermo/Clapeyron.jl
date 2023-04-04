
function rachfordrice_β0(K,z,β0 = nothing)
    g0 = dot(z, K) - 1.
    g1 = 1. - sum(zi/Ki for (zi,Ki) in zip(z,K))
    singlephase = false
    _1 = one(g1)
    _0 = zero(g1)
    if g0 < 0
        β = _0
        singlephase = true
    elseif g1 > 0
        β = _1
        singlephase = true
    end

    β0 !== nothing && return β0,singlephase
    βmin =  _0
    βmax = _1
    for i in eachindex(K)
        Ki,zi = K[i],z[i]
        if Ki > 1
            βmin = min(βmin,(Ki*zi - 1)/(Ki - 1))
        end

        if Ki < 1
            βmax = max(βmax,(1 - zi)/(1 - Ki))
        end
    end
    
    β = (βmax + βmin)/2

    return β,singlephase
end

function rachfordrice(K, z; β0=nothing, non_inx=FillArrays.Fill(false,length(z)), non_iny=non_inx)
    # Function to solve Rachdord-Rice mass balance
    β,singlephase = rachfordrice_β0(K,z,β0)
    
    if length(z) <= 4 && all(Base.Fix2(>,0),z) && all(!,non_inx) && all(!,non_iny) && !singlephase
        return rr_vle_vapor_fraction_exact(K,z)
    end

    #halley refinement
    if !singlephase
        return rr_flash_refine(K,z,β,non_inx,non_iny)
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
             [idiag] += 1/ny[i]
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

    ϕx = log.(x) .+ lnϕx
    ϕy = log.(y) .+ lnϕy

    # to avoid NaN in Gibbs energy
    for i in eachindex(z)
        non_iny[i] && (ϕy[i] = 0.)
        non_inx[i] && (ϕx[i] = 0.)
    end

    if G !== nothing
        # Computing Gibbs Energy gradient
        if !all_equilibria
            G .= (ϕy .- ϕx)[in_equilibria]
        else
            G .= ϕy .- ϕx
        end
    end

    if F !== nothing
        # Computing Gibbs Energy
        FO = dot(ny,ϕy) + dot(nx,ϕx)
        return FO
    end
end

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
- `noncondensables` = arrays with names (strings) of components non allowed on the liquid phase. Not allowed with `lle` equilibria
- `nonvolatiles` = arrays with names (strings) of components non allowed on the vapour phase. Not allowed with `lle` equilibria

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

function index_reduction(m::MichelsenTPFlash,idx::AbstractVector)
    equilibrium,K0,x0,y0,v0,K_tol,ss_iters,second_order,noncondensables,nonvolatiles = m.equilibrium,m.K0,m.x0,m.y0,m.v0,m.K_tol,m.ss_iters,m.second_order,m.noncondensables,m.nonvolatiles
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return MichelsenTPFlash(;equilibrium,K0,x0,y0,v0,K_tol,ss_iters,second_order,noncondensables,nonvolatiles)
end

numphases(::MichelsenTPFlash) = 2

function MichelsenTPFlash(;equilibrium = :vle,
                        K0 = nothing, 
                        x0 = nothing,
                        y0=nothing,
                        v0=nothing,
                        K_tol = sqrt(eps(Float64)),
                        ss_iters = 21,
                        nacc = 5,
                        second_order = false,
                        noncondensables = nothing,
                        nonvolatiles = nothing)
    !(is_vle(equilibrium) | is_lle(equilibrium)) && throw(error("invalid equilibrium specification for MichelsenTPFlash"))
    if K0 == x0 == y0 === v0 == nothing #nothing specified
        is_lle(equilibrium) && throw(error("""
        You need to provide either an initial guess for the partion constant K
        or for compositions of x and y for LLE"""))
        T = nothing
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
    if is_lle(equilibrium)
        if !isnothing(nonvolatiles) && length(nonvolatiles) > 0
            throw(error("LLE equilibria does not support setting nonvolatiles"))
        end
    
        if !isnothing(noncondensables) && length(noncondensables) > 0
            throw(error("LLE equilibria does not support setting noncondensables"))
        end
    end

    return MichelsenTPFlash{T}(equilibrium,K0,x0,y0,v0,K_tol,ss_iters,nacc,second_order,noncondensables,nonvolatiles)
end

is_vle(method::MichelsenTPFlash) = is_vle(method.equilibrium)
is_lle(method::MichelsenTPFlash) = is_lle(method.equilibrium)

function tp_flash_impl(model::EoSModel,p,T,z,method::MichelsenTPFlash)
    x,y,β =  tp_flash_michelsen(model,p,T,z;equilibrium = method.equilibrium, K0 = method.K0,
            x0 = method.x0, y0 = method.y0, vol0 = method.v0,
            K_tol = method.K_tol,itss = method.ss_iters, nacc=method.nacc,
            second_order = method.second_order,
            non_inx_list=method.noncondensables, non_iny_list=method.nonvolatiles, 
            reduced = true)
    
    G = (gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/R̄/T

    X = hcat(x,y)'
    nvals = X.*[1-β
                β] .* sum(z)
    return (X, nvals, G)
end

#updates lnK, returns lnK,volx,voly, gibbs if β != nothing
function update_K!(lnK,model,p,T,x,y,volx,voly,phasex,phasey,β = nothing)
    lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
    lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
    lnK .= lnϕx - lnϕy
    gibbs = zero(eltype(lnK))
    if β !== nothing
        for i in 1:length(y)
            if iny[i]
                gibbs += β*y[i]*log(y[i] + lnϕy[i])
            end
            if inx[i]
                gibbs += (1-β)*x[i]*log(x[i] + lnϕx[i])
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

function tp_flash_michelsen(model::EoSModel, p, T, z; equilibrium=:vle, K0=nothing,
                                     x0=nothing, y0=nothing, vol0=(nothing, nothing),
                                     K_tol=1e-8, itss=21, nacc=5, second_order=false,
                                     non_inx_list=String[], non_iny_list=String[], reduced=false)


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
    if !isnothing(K0)
        K = 1. * K0
        lnK = log.(K)
    elseif !isnothing(x0) && !isnothing(y0)
        x = x0
        y = y0
        lnK,volx,voly,_ = update_K!(lnK,model,p,T,x,y,volx,voly,phasex,phasey,nothing)
        K = exp.(lnK)
    elseif is_vle(equilibrium)
        # Wilson Correlation for K
        K = wilson_k_values(model,p,T)
        lnK = log.(K)
    else
        err() = @error("""You need to provide either an initial guess for the partion constant K
                        or for compositions of x and y for LLE""")
        err()
    end

    _1 = one(p+T+first(z))
    # Initial guess for phase split
    β,singlephase = rachfordrice_β0(K0,z)
    #=TODO:
    there is a method used in TREND that tries to obtain adequate values of K
    in the case of incorrect initialization.
    =#

    # Stage 1: Successive Substitution
    error_lnK = _1
    it = 0

    x = similar(z)
    y = similar(z)
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

    while error_lnK > K_tol && it < itss && !singlephase
        it += 1
        lnK_old = lnK .* _1

        β = rachfordrice(K, z; β0=β, non_inx=non_inx, non_iny=non_iny)

        singlephase = !(0 < β < 1) #rachford rice returns 0 or 1 if it is single phase.
        x,y = update_rr!(K,β,z,x,y,non_inx,non_iny)
        # Updating K's
        lnK,volx,voly,gibbs = update_K!(lnK,model,p,T,x,y,volx,voly,phasex,phasey,β)
    
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
            lnK_dem,volx_dem,voly_dem,gibbs_dem = update_K!(lnK_dem,model,p,T,x_dem,y_dem,volx_dem,voly_dem,phasex,phasey,β)
            # only accelerate if the gibbs free energy is reduced
            if gibbs_dem < gibbs 
                lnK .= _1 * lnK_dem
                volx = _1 * volx_dem
                voly = _1 * voly_dem
                β = _1 * β_dem
            end
        end

        K .= exp.(lnK)

        # Computing error
        # error_lnK = sum((lnK .- lnK_old).^2)
        error_lnK = dnorm(lnK,lnK_old,1)
    end

    # Stage 2: Minimization of Gibbs Free Energy
    vcache = Ref((volx, voly))

    if error_lnK > K_tol && it == itss &&  !singlephase
        # println("Second order minimization")
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

    if singlephase
        β = zero(β)/zero(β)
        x .= z
        y .= z
    end

    if !reduced
        x = index_expansion(x,z_nonzero)
        y = index_expansion(y,z_nonzero)
    end
    vx,vy = vcache[]
    if vx < vy #sort by increasing volume
        return x, y, β
    else
        return y, x, 1 - β
    end
end

export MichelsenTPFlash
