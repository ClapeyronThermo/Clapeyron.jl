#=
function rachfordrice(β, K, z)
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

    # Solving the phase fraction β using Halley's method
    it = 0
    error = 1.
    while error > 1e-8 && it < 10 &&  ~singlephase
        it = it + 1
        D = 1. .+ β.*K1
        KD = K1./D
        FO = dot(z, KD)
        dFO = - dot(z, KD.^2)
        d2FO = 2. *dot(z, KD.^3)
        dβ = - (2*FO*dFO)/(2*dFO^2-FO*d2FO)
        β = β + dβ
        error = abs(dβ)
    end
    return β, D, singlephase
end

=#
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
        ∂ϕy ./= nxsum
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
- `second_order` = wheter to solve the gibbs energy minimization using the analytical hessian or not
"""
struct MichelsenTPFlash{T} <: TPFlashMethod
    equilibrium::Symbol
    K0::Union{Vector{T},Nothing}
    x0::Union{Vector{T},Nothing}
    y0::Union{Vector{T},Nothing}
    v0::Union{Tuple{T,T},Nothing}
    K_tol::Float64
    ss_iters::Int
    second_order::Bool
end

function index_reduction(m::MichelsenTPFlash,idx::AbstractVector)
    equilibrium,K0,x0,y0,v0,K_tol,ss_iters,second_order = m.equilibrium,m.K0,m.x0,m.y0,m.v0,m.K_tol,m.ss_iters,m.second_order
    K0 !== nothing && (K0 = K0[idx])
    x0 !== nothing && (x0 = x0[idx])
    y0 !== nothing && (y0 = y0[idx])
    return MichelsenTPFlash(;equilibrium,K0,x0,y0,v0,K_tol,ss_iters,second_order)
end

numphases(::MichelsenTPFlash) = 2

function MichelsenTPFlash(;equilibrium = :vle,K0 = nothing, x0 = nothing,y0=nothing,v0=nothing,K_tol = eps(Float64),ss_iters = 10,second_order = false)
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
    return MichelsenTPFlash{T}(equilibrium,K0,x0,y0,v0,K_tol,ss_iters,second_order)
end

is_vle(method::MichelsenTPFlash) = is_vle(method.equilibrium)
is_lle(method::MichelsenTPFlash) = is_lle(method.equilibrium)

function tp_flash_impl(model::EoSModel,p,T,z,method::MichelsenTPFlash)
    x,y,β =  tp_flash_michelsen(model,p,T,z;equilibrium = method.equilibrium, K0 = method.K0,
                        x0 = method.x0, y0 = method.y0, vol0 = method.v0,
                        K_tol = method.K_tol,itss = method.ss_iters,second_order = method.second_order,
                        reduced = true)


    G = (gibbs_free_energy(model,p,T,x)*(1-β)+gibbs_free_energy(model,p,T,y)*β)/R̄/T

    X = hcat(x,y)'
    nvals = X.*[1-β
                β] .* sum(z)
    return (X, nvals, G)
end

function tp_flash_michelsen(model::EoSModel, p, T, z; equilibrium=:vle, K0=nothing,
                            x0=nothing, y0=nothing, vol0=(nothing, nothing),
                            K_tol=1e-16, itss=10, second_order=false, reduced = false)
    # Function to compute two phase flash at given temperature, pressure and
    # global composition
    # p = Pressure
    # T = Temperature
    # z = global composition array
    # equilibrium = equilibrium type ":vle" for liquid vapor equilibria, ":lle" for liquid liquid equilibria
    # K0 = optional, initial guess for the constants K
    # x0 = optional, initial guess for the composition of phase x
    # y0 = optional, initial guess for the composition of phase y
    # vol0 = optional, initial guesses for phase x and phase y volumes
    # K_tol = tolerance to stop the calculation
    # itss = number of Successive Substitution iterations to perform
    # second_order = wheter to solve the gibbs energy minimization using the analycal hessian or not
    #reduced = if the model has been striped of nonzero values
    # out = phase x composition, phase y composition, phase split fraction
    #reduce model
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
    _1 = one(p+T+first(z))
    # Initial guess for phase split
    βmin = max(0., minimum(((K.*z .- 1) ./ (K .-  1.))[K .> 1]))
    βmax = min(1., maximum(((1 .- z) ./ (1. .- K))[K .< 1]))
    β = _1*(βmin + βmax)/2
    # Stage 1: Successive Substitution
    it = 0
    error = _1
    singlephase = false
    lnK_old = copy(K) .* _1
    x = similar(z)
    y = similar(z)
    while error > K_tol && it < itss
        it += 1
        lnK_old .= lnK
        # Solving Rachford-Rice Eq.
        β = rr_vle_vapor_fraction(K, z)
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
        K .= exp.(lnK)
        # Computing error
        error = dnorm(lnK,lnK_old,1)
    end
    vt = volx,voly
    vcache = Ref{typeof(vt)}(vt)


    # Stage 2: Minimization of Gibbs Free Energy
    if error > K_tol && it == itss &&  ~singlephase
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
        fill!(x,z)
        fill!(y,z)
    end

    if !reduced
        x = index_expansion(x,z_nonzero)
        y = index_expansion(y,z_nonzero)
    end

    return x, y, β
end

export MichelsenTPFlash