function rachfordrice(K, z; β0=nothing)
    # Function to solve Rachdord-Rice mass balance
    K1 = K .- 1.
    g0 = dot(z, K) - 1.
    g1 = 1. - dot(z, 1. ./ K)
    singlephase = false

    # Checking if the given K and z have solution
    if g0 < 0
        β = 0.
        D = fill!(similar(z), 1)
        singlephase = true
    elseif g1 > 0
        β = 1.
        D = 1 .+ K1
        singlephase = true
    end

    βmin = max(0., minimum(((K.*z .- 1) ./ (K .-  1.))[K .> 1]))
    βmax = min(1., maximum(((1 .- z) ./ (1. .- K))[K .< 1]))

    if isnothing(β0)
        β = (βmax + βmin)/2
    else
        β = 1. * β0
    end

    # Solving the phase fraction β using Halley's method
    it = 0
    error_β = 1.
    error_FO = 1.
    while error_β > 1e-8 && error_FO > 1e-8 && it < 10 &&  ~singlephase
        it = it + 1
        D = 1. .+ β.*K1
        KD = K1./D
        FO = dot(z, KD)
        dFO = - dot(z, KD.^2)
        d2FO = 2. * dot(z, KD.^3)
        dβ = - (2*FO*dFO)/(2*dFO^2-FO*d2FO)

        # restricted β space
        if FO < 0.
            βmax = β
        elseif FO > 0.
            βmin = β
        end

        #updatind β
        βnew =  β + dβ
        if βmin < βnew && βnew < βmax
            β = βnew
        else
            dβ = (βmin + βmax) / 2 - β
            β = dβ + β
        end

        error_β = abs(dβ)
        error_FO = abs(FO)

    end
    return β
end


function gibbs_obj!(model::EoSModel, p, T, z, phasex, phasey, ny, vcache; F=nothing, G=nothing)
    # Objetive Function to minimize the Gibbs Free Energy
    # It computes the Gibbs free energy and its gradient
    nx = z .- ny
    x = nx ./ sum(nx)
    y = ny ./ sum(ny)

    # Volumes are set from local cache to reuse their values for following
    # Iterations
    volx,voly = vcache[]

    lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
    lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)

    ϕx = log.(x) .+ lnϕx
    ϕy = log.(y) .+ lnϕy

    #volumes are stored in the local cache
    vcache[] = (volx,voly)
    if G !== nothing
        # Computing Gibbs Energy gradient
        G .= ϕy .- ϕx
    end

    if F !== nothing
        # Computing Gibbs Energy
        FO = dot(ny,ϕy) + dot(nx,ϕx)
        return FO
    end

end


function dgibbs_obj!(model::EoSModel, p, T, z, phasex, phasey, ny, vcache; F=nothing, G=nothing, H=nothing)
    # Objetive Function to minimize the Gibbs Free Energy
    # It computes the Gibbs free energy, its gradient and its hessian
    nx = z .- ny
    nxsum = sum(nx)
    nysum = sum(ny)
    x = nx ./ nxsum
    y = ny ./ nysum

    # Volumes are set from local cache to reuse their values for following
    # Iterations
    volx,voly = vcache[]

    if H !== nothing
        # Computing Gibbs Energy Hessian
        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y; phase=phasey, vol0=voly)

        ∂ϕx,∂ϕy = ∂lnϕ∂nx,∂lnϕ∂ny

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
    ϕx = log.(x) .+ lnϕx
    ϕy = log.(y) .+ lnϕy

    if G !== nothing
        # Computing Gibbs Energy gradient
        G .= ϕy .- ϕx
    end

    if F != nothing
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
    modified = false #if there is any nonvolatiles/noncondensables  
    if !isnothing(method.nonvolatiles) && length(method.nonvolatiles) > 0
        modified = true
    end

    if !isnothing(method.noncondensables) && length(method.noncondensables) > 0
        modified = true
    end


    if !modified
        x,y,β =  tp_flash_michelsen(model,p,T,z;equilibrium = method.equilibrium, K0 = method.K0,
                            x0 = method.x0, y0 = method.y0, vol0 = method.v0,
                            K_tol = method.K_tol,itss = method.ss_iters, nacc=method.nacc,
                            second_order = method.second_order,
                            reduced = true)
    else
        x,y,β =  tp_flash_michelsen_modified(model,p,T,z;equilibrium = method.equilibrium, K0 = method.K0,
                        x0 = method.x0, y0 = method.y0, vol0 = method.v0,
                        K_tol = method.K_tol,itss = method.ss_iters, nacc=method.nacc,
                        second_order = method.second_order,
                        non_inx_list=method.noncondensables, non_iny_list=method.nonvolatiles, reduced = true)
    end


    G = (gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/R̄/T

    X = hcat(x,y)'
    nvals = X.*[1-β
                β] .* sum(z)
    return (X, nvals, G)
end

function tp_flash_michelsen(model::EoSModel, p, T, z; equilibrium=:vle, K0=nothing,
                            x0=nothing, y0=nothing, vol0=(nothing, nothing),
                            K_tol=1e-8, itss=21, nacc=5, second_order=false, reduced = false)

    if !reduced
        model_full,z_full = model,z
        model, z_nonzero = index_reduction(model_full,z_full)
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

    # Computing the initial guess for the K vector
    if ~isnothing(K0)
        K = 1. * K0
        lnK = log.(K)
    elseif ~isnothing(x0) && ~isnothing(y0)
        x = x0
        y = y0
        lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
        lnK = lnϕx - lnϕy
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

    # type stability
    _1 = one(p+T+first(z))

    # Initial guess for phase split
    β0 = nothing
    β = 0.5 * _1 # this is just to have β as a global variable

    # Stage 1: Accelerated Successive Substitution
    it = 0
    error_lnK = copy(_1)
    singlephase = false
    lnK_old = lnK .* _1
    x = similar(z)
    y = similar(z)

    itacc = 0
    lnK3 = similar(lnK)
    lnK4 = similar(lnK)
    lnK5 = similar(lnK)
    lnK_dem = similar(lnK)

    gibbs = copy(_1)
    gibbs_dem = copy(_1)

    while error_lnK > K_tol && it < itss
        itacc += 1
        it += 1
        lnK_old .= lnK

        # solving RR's equations using Halley's method
        β = rachfordrice(K, z; β0=β0)

        singlephase = !(0 <= β <= 1)
        # Recomputing phase composition
        x = rr_flash_liquid!(x,K,z,β)
        y .= K .* x
        x ./= sum(x)
        y ./= sum(y)
        # Updating K's
        lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
        lnK .= lnϕx .- lnϕy

        # computing current Gibbs free energy 
        gibbs = β * sum(y .* (log.(y) .+ lnϕy)) 
        gibbs += (1. - β) * sum(x .* (log.(x) .+ lnϕx))
        
        # acceleration step
        if itacc == (nacc - 2)
            lnK3 = 1. * lnK
        elseif itacc == (nacc - 1)
            lnK4 = 1. * lnK
        elseif itacc == nacc
            itacc = 0
            lnK5 = 1. * lnK
            # acceleration using DEM (1 eigenvalues)
            lnK_dem = dem(lnK5, lnK4, lnK3)
            K_dem = exp.(lnK_dem)
            # requires to solve the RR problem again (to compute the 
            # extrapolated Gibbs free energy)
            β_dem = rachfordrice(K_dem, z; β0=β)
            x_dem = similar(x)
            x_dem = rr_flash_liquid!(x_dem,K_dem,z,β_dem)
            y_dem = x_dem .* K_dem
            x_dem ./= sum(x_dem)
            y_dem ./= sum(y_dem)

            lnϕx_dem, volx_dem = lnϕ(model, p, T, x_dem; phase=phasex, vol0=volx)
            lnϕy_dem, voly_dem = lnϕ(model, p, T, y_dem; phase=phasey, vol0=voly)

            # computing the extrapolated Gibbs free energy
            gibbs_dem = β_dem * sum(y_dem .* (log.(y_dem) .+ lnϕy_dem)) 
            gibbs_dem += (1. - β_dem) * sum(x_dem .* (log.(x_dem) .+ lnϕx_dem))

            # only accelerate if the gibbs free energy is reduced
            if gibbs_dem < gibbs 
                lnK = _1 * lnK_dem
                volx = _1 * volx_dem
                voly = _1 * voly_dem
                β = _1 * β_dem
            end
        end
        
        K .= exp.(lnK)
        # updating future initial guess for the phase fraction 
        β0 = _1 * β
        # Computing error
        # error_lnK = dnorm(lnK,lnK_old,1)
        error_lnK = sum(abs.(lnK - lnK_old))
    end

    vt = volx,voly
    vcache = Ref{typeof(vt)}(vt)

    
    # Stage 2: Minimization of Gibbs Free Energy
    if error_lnK > K_tol && it == itss &&  ~singlephase
        ny = β*y
        # minimizing Gibbs Free Energy
        if second_order
            dfgibbs!(F, G, H, ny) = dgibbs_obj!(model, p, T, z, phasex, phasey,
                                             ny,vcache; F=F, G=G, H=H)
            #sol = Optim.optimize(only_fgh!(dfgibbs!), ny, Optim.Newton())
            sol = Solvers.optimize(Solvers.only_fgh!(dfgibbs!), ny, LineSearch(Newton()))
        else
            fgibbs!(F, G, ny) = gibbs_obj!(model, p, T, z, phasex, phasey,
                                           ny,vcache; F, G=G)
            sol = Solvers.optimize(Solvers.only_fg!(fgibbs!), ny, LineSearch(BFGS()))
        end
        # Converting from moles to mole fractions
        # ny = sol.minimizer
        ny = Solvers.x_sol(sol)
        β = sum(ny)
        nx = z .- ny
        x = nx ./ sum(nx)
        y = ny ./ β

    end

    if singlephase
        β = zero(β)/zero(β)
        # Gustavo: the fill! function was giving an error
        # fill!(x,z)
        # fill!(y,z)
        x .= z
        y .= z
    end

    if !reduced
        x = index_expansion(x,z_nonzero)
        y = index_expansion(y,z_nonzero)
    end

    return x, y, β
end

export MichelsenTPFlash
