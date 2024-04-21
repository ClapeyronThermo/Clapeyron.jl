"""
    GuptaMultiphaseFlash{T, T2}(;K0 = nothing,β0 = nothing,rtol= 1e-10,atol = 1e-10,max_iters = 100)
The optimizer will stop at `max_iters` evaluations, when the absolute tolerance is less than `atol` or when the relative tolerance is less than `rtol`
"""
Base.@kwdef struct GuptaMultiphaseFlash{T} <: TPFlashMethod  
    K0::T = nothing
    β0 = nothing
    rtol::Float64 = 1e-10
    atol::Float64 = 1e-10
    max_iters::Int = 100
    spec::Symbol = :unknown
end


index_reduction(flash::GuptaMultiphaseFlash{Nothing},idx::AbstractVector) = flash

function index_reduction(flash::GuptaMultiphaseFlash,idx::AbstractVector)
    K02 = flash.K0[idx]
    return GuptaMultiphaseFlash(K02,flash.β0,flash.rtol,flash.atol,flash.max_iters,flash.spec)
end

function tp_flash_impl(model::EoSModel, p, T, n, method::GuptaMultiphaseFlash)
    z = normalize(n, 1)
    nc = length(model)
    #TODO make type-stable. Needs some semantics for handling β0, K0
    # T = eltype(promote(method.β0, method.K0))

    # Initialise arrays
    if method.β0 isa Nothing
        β = normalize(ones(Float64, 2), 1)
        np = 2
    else
        β = method.β0
        np = length(β0)
    end


    if method.K0 isa Nothing
        @warn "If K0 is not specified, initialisation with more than 2 phases is not supported"
        Kᵂ = wilson_k_values(model, p, T)
        K = ones(eltype(Kᵂ), np, nc)
        K[1, :] .= eltype(Kᵂ)(1)
        K[2, :] .= Kᵂ
    else
        K = method.K0
    end

    # Iteration vector
    x = deepcopy(β[2:np])
    for i = 2:np
        if x[i-1] == 0
            x[i-1] = -1
        end
    end

    iter = 0
    converge_norm = 100.0
    x_mat = zeros(np, nc)
    ϕ_mat = zeros(np, nc)
    while converge_norm > method.atol && iter < method.max_iters
        iter += 1
        x_old = deepcopy(x)

        # Calculate mole fractions from K-values
        update_x_mat!(x_mat, x, K, z, np, nc)

        # Calculate fugacity coefficients
        for (j, x_row) in enumerate(eachrow(x_mat))
            ϕ_mat[j, :] .= fugacity_coefficient(model, p, T, x_row)
        end

        # Calculate K-values
        for j = 1:np
            for i = 1:nc
                K[j, i] = ϕ_mat[1, i]/ϕ_mat[j, i]
            end
        end

        # F = zeros(np-1)
        # Non-linear root-finding problem in Rⁿᵖ⁻¹
        obj_func!(F, x) = Gupta_obj_func!(F, x, K, z, np, nc)

        # Chunk size of 1 as for 2 phase problems this drops to a 1D problem. Ideally should be able to offload to ForwardDiff inbuilt heuristic??
        res = Solvers.nlsolve(obj_func!, x, TrustRegion(Newton(), NWI()), NEqOptions(), ForwardDiff.Chunk{1}())
        x = Solvers.x_sol(res)

        # Transform to phase fraction and stability values
        β, θ = unpack_xvec(x, np)
        for i = 2:np
            if θ[i] > 0
                x[i-1] = -θ[i]
            else
                x[i-1] = β[i] 
            end
        end

        # Calculate convergence norm (orig. code used 1-norm)
        converge_norm = norm(x .- x_old, Inf)
    end

    iter == method.max_iters && @warn("Did not converge in $max_iters iterations, with a residual of $converge_norm")
    #TODO Test if shift in phase order required

    @bp
    return x_mat, β
end

function Gupta_obj_func!(F, x, K, z, np, nc)
    β, θ = unpack_xvec(x, np)
    T = eltype(x)

    # Calculate objective function
    # F = zeros(T, np-1)
    for j = 2:np
        sumat = T(0)
        for i = 1:nc
            sum_int = T(0)
            for k = 2:np
                sum_int = sum_int + β[k]*(K[k, i]*exp(θ[k]) - 1)
            end
            sumat = sumat + z[i]*(K[j, i]*exp(θ[j])-1)/(1+sum_int)
        end
        F[j-1] = sumat
    end

    return F
end

function unpack_xvec(x, np)
    T1 = eltype(x)
    # Variable vectors
    _0 = T1(0)
    _1 = T1(1)
    ϵ = eps(T1)

    β = zeros(T1, np)
    θ = zeros(T1, np)

    for i = 2:np
        if x[i-1] < 0
            θ[i] = -x[i-1]
        else
            β[i] = x[i-1]
        end
    end

    # Take first phase as reference, requires β > 0 & θ = 0
    θ[1] = _0
    β[1] = _1 - sum(β)

    # Round values < ϵ
    #? Could this be moved to the initial read-in of x?
    β[β .< ϵ] .= _0
    θ[θ .< ϵ] .= _0

    # Normalise β
    β = normalize(abs.(β), 1)
    return β, θ
end

function update_x_mat!(x_mat, x, K, z, np, nc)
    β, θ = unpack_xvec(x, np)

    for i = 1:nc
        sumat = 0.0
        for j = 2:np
            sumat = sumat + (K[j, i]*exp(θ[j]) - 1)*β[j]
        end
        x_mat[1, i] = z[i]/(1+sumat)
    end

    for i = 1:nc
        for j = 2:np
            x_mat[j, i] = x_mat[1,i]*K[j,i]*exp(θ[j])
        end
    end

    return x_mat
end

export GuptaMultiphaseFlash