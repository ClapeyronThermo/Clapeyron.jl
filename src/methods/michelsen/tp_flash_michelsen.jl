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

function gibbs_obj!(model::EoSModel, p, T, z, z_notzero, phasex, phasey, ny; F=nothing, G=nothing)
    # Objetive Function to minimize the Gibbs Free Energy
    # It computes the Gibbs free energy and its gradient
    x = fill!(similar(z), 0)
    y = fill!(similar(z), 0)

    nx = z[z_notzero] .- ny
    x[z_notzero] = nx ./ sum(nx)
    y[z_notzero] = ny ./ sum(ny)

    # Volumes are set to global variables to reuse their values for following
    # Iterations
    global volx
    global voly

    lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
    lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)

    ϕx = log.(x[z_notzero]) + lnϕx[z_notzero]
    ϕy = log.(y[z_notzero]) + lnϕy[z_notzero]

    if G != nothing
        # Computing Gibbs Energy gradient
        G[:] = ϕy - ϕx
    end

    if F != nothing
        # Computing Gibbs Energy
        FO = sum(ny.*ϕy[z_notzero] + nx.*ϕx[z_notzero])
        return FO
    end

end


function dgibbs_obj!(model::EoSModel, p, T, z, z_notzero, phasex, phasey, ny; F=nothing, G=nothing, H=nothing)
    # Objetive Function to minimize the Gibbs Free Energy
    # It computes the Gibbs free energy, its gradient and its hessian
    x = fill!(similar(z), 0)
    y = fill!(similar(z), 0)

    ncomponents = length(ny)
    nx = z[z_notzero] .- ny
    nxsum = sum(nx)
    nysum = sum(ny)
    x[z_notzero] = nx ./ nxsum
    y[z_notzero] = ny ./ nysum

    # Volumes are set to global variables to reuse their values for following
    # Iterations
    global volx
    global voly

    if H != nothing
        # Computing Gibbs Energy Hessian
        lnϕx, ∂lnϕ∂nx, ∂lnϕ∂Px, volx = ∂lnϕ∂n∂P(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, ∂lnϕ∂ny, ∂lnϕ∂Py, voly = ∂lnϕ∂n∂P(model, p, T, y; phase=phasey, vol0=voly)

        eye = Matrix{Float64}(Identity, ncomponents, ncomponents)
        ∂ϕx = eye./nx .- 1/nxsum .+ ∂lnϕ∂nx[z_notzero, z_notzero]/nxsum
        ∂ϕy = eye./ny .- 1/nysum .+ ∂lnϕ∂ny[z_notzero, z_notzero]/nysum
        H[:, :] = ∂ϕx + ∂ϕy
    else
        lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
    end

    ϕx = log.(x[z_notzero]) + lnϕx[z_notzero]
    ϕy = log.(y[z_notzero]) + lnϕy[z_notzero]

    if G != nothing
        # Computing Gibbs Energy gradient
        G[:] = ϕy - ϕx
    end

    if F != nothing
        # Computing Gibbs Energy
        FO = sum(ny.*ϕy[z_notzero] + nx.*ϕx[z_notzero])
        return FO
    end

end

"""
    MichelsenTPFlash{T}(;K0 = nothing,rtol= 1e-10,atol = 1e-10,max_iters = 100)

Method to solve non-reactive multicomponent flash problem by Michelsen's method.

Only two phases are supported. if `K0` is `nothing`, it will be calculated via the Wilson correlation.

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
    K_tol::Float64 = eps(Float64)
    ss_iters::Int = 10
    second_order::Bool = false
end

function MichelsenTPFlash(;equilibrium = :vle,K0 = nothing, x0 = nothing,y0=nothing,v0=nothing,K_tol = eps(Float64),ss_iters = 10,second_order = false)
    !(is_vle(equilibrium) | is_lle(equilibrium)) && throw(error("invalid equilibrium specification for MichelsenTPFlash"))
    if K0 == x0 == y0 === v0 == nothing #nothing specified
        equilibrium == :lle && throw(error("""
        You need to provide either an initial guess for the partion constant K
        or for compositions of x and y for LLE"""))   
        T = Float64
    else
        if !isnothing(K0) & isnothing(x0) & isnothing(y0) #K0 specified
            T = eltype(K0)
        elseif isnothing(K0) & isnothing(x0) & isnothing(y0)  #x0, y0 specified
            T = eltype(x0)
        else
            throw(error("invalid specification of initial points"))
        end
     end
    return MichelsenTPFlash{T}(equilibrium,K0,x0,y0,v0,K_tol,ss_iters,second_order)
end

function tp_flash_michelsen(model::EoSModel, p, T, z; equilibrium=:vle, K0=nothing,
                            x0=nothing, y0=nothing, vol0=[nothing, nothing],
                            K_tol=1e-16, itss=10, second_order=false)
    # Function to compute two phase flash at given temperature, pressure and
    # global composition
    # p = Pressure
    # T = Temperature
    # z = global composition array
    # equilibrium = equilibrium type ":lv" for liquid vapor equilibria, ":lle" for liquid liquid equilibria
    # K0 = optional, initial guess for the constants K
    # x0 = optional, initial guess for the composition of phase x
    # y0 = optional, initial guess for the composition of phase y
    # vol0 = optional, initial guesses for phase x and phase y volumes
    # K_tol = tolerance to stop the calculation
    # itss = number of Successive Substitution iterations to perform
    # second_order = wheter to solve the gibbs energy minimization using the analycal hessian or not

    # out = phase x composition, phase y composition, phase split fraction


    if is_vle(equilibrium)
        phasex = :liquid
        phasey = :vapor
    elseif is_lle(equilibrium)
        phasex = :liquid
        phasey = :liquid
    end

    # Setting the initial guesses for volumes
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
        # Check this function, it didnt work with SAFT-γ-Mie
        K = wilson_k_values(model,p,T)
    else
        err() = @error("""You need to provide either an initial guess for the partion constant K
                        or for compositions of x and y for LLE""")
        err()
    end

    # Initial guess for phase split
    βmin = max(0., minimum(((K.*z .- 1) ./ (K .-  1.))[K .> 1]))
    βmax = min(1., maximum(((1 .- z) ./ (1. .- K))[K .< 1]))
    β = (βmin + βmax)/2

    # Stage 1: Successive Substitution
    it = 0
    error = 1.
    singlephase = false
    while error > K_tol && it < itss
        it += 1
        lnK_old = lnK
        # Solving Rachford-Rice Eq.
        β, D, singlephase = rachfordrice(β, K, z)
        # Recomputing phase composition
        x = z ./ D
        y = x .* K
        x ./= sum(x)
        y ./= sum(y)
        # Updating K's
        lnϕx, volx = lnϕ(model, p, T, x; phase=phasex, vol0=volx)
        lnϕy, voly = lnϕ(model, p, T, y; phase=phasey, vol0=voly)
        lnK = lnϕx - lnϕy
        K = exp.(lnK)
        # Computing error
        error = sum(abs.(lnK - lnK_old))
    end

    # Stage 2: Minimization of Gibbs Free Energy
    if error > K_tol && it == itss &&  ~singlephase
        global volx
        global voly
        z_notzero = z .> 0.
        ny = β*y[z_notzero]
        # minimizing Gibbs Free Energy
        if second_order
            dfgibbs!(F, G, H, ny) = dgibbs_obj!(model, p, T, z, z_notzero, phasex, phasey,
                                             ny; F=F, G=G, H=H)
            sol = Optim.optimize(only_fgh!(dfgibbs!), ny, Optim.Newton())
        else
            fgibbs!(F, G, ny) = gibbs_obj!(model, p, T, z, z_notzero, phasex, phasey,
                                           ny; F, G=G)
            sol = Optim.optimize(only_fg!(fgibbs!), ny, Optim.BFGS())
        end
        # Converting from moles to mole fractions
        ny = sol.minimizer
        x = fill!(similar(z), 0)
        y = fill!(similar(z), 0)
        β = sum(ny)
        nx = z[z_notzero] .- ny
        x[z_notzero] = nx ./ sum(nx)
        y[z_notzero] = ny ./ β
    end
    return x, y, β
end
