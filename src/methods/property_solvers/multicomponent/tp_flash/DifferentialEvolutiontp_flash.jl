#=
Original code by Thomas Moore
@denbigh
included in https://github.com/ClapeyronThermo/Clapeyron.jl/pull/56
=#
"""
    DETPFlash(; numphases = 2,
    max_steps = 3000(numphases-1),
    population_size = 50,
    time_limit = Inf,
    seed = 1,
    stagnation_evals = 0,
    stagnation_tol = 0.0,
    backend = :sass,
    verbose = false,
    logspace = false,
    equilibrium = :auto)

Method to solve a non-reactive multicomponent TP-flash problem by global optimization of the Gibbs free energy.

This implementation uses an ask–tell metaheuristic backend. The current default backend is `:sass`
(Self-adaptive spherical search, SASS). The optimizer stops when it reaches `max_steps` function
evaluations, exceeds `time_limit` seconds, or triggers stagnation-based early stopping controlled by
`stagnation_evals` and `stagnation_tol`.

Note: the file name is historical; the original implementation used differential evolution.

User must assume a number of phases, `numphases`. If true number of phases is smaller than numphases, model should predict either (a) identical composition in two or more phases, or (b) one phase with negligible total number of moles. If true number of phases is larger than numphases, a thermodynamically unstable solution will be predicted.

The `equilibrium` keyword allows to restrict the search of phases to just liquid-liquid equilibria (`equilibrium = :lle`). the default searches for liquid and gas phases.
"""
@kwdef struct DETPFlash <: TPFlashMethod
    numphases::Int = 2
    max_steps::Int = 3000(numphases-1)
    population_size::Int = 50
    time_limit::Float64 = Inf
    seed::Int = 1
    stagnation_evals::Int = 0
    stagnation_tol::Float64 = 0.0
    backend::Symbol = :sass
    verbose::Bool = false
    logspace::Bool = false
    equilibrium::Symbol = :auto
end

index_reduction(flash::DETPFlash,z) = flash

#z is the original feed composition, x is a matrix with molar fractions, n is a matrix with molar amounts
function partition!(dividers,n,x,nvals)
    numphases, numspecies = size(x)
    @inbounds for j = 1:numspecies
        nj = n[j]
        remaining = nj
        for i = 1:(numphases - 1)
            nij = dividers[i, j]*remaining
            nvals[i, j] = nij
            remaining -= nij
        end
        nvals[numphases, j] = remaining
    end
    #Calculate mole fractions xij
    @inbounds for i = 1:numphases
        s = zero(eltype(nvals))
        @simd for j = 1:numspecies
            s += nvals[i, j]
        end
        if iszero(s)
            fill!(@view(x[i, :]), zero(eltype(x)))
        else
            invs = one(s)/s
            @simd for j = 1:numspecies
                x[i, j] = nvals[i, j]*invs
            end
        end
    end
end

function tp_flash_impl(model::EoSModel, p, T, n, method::DETPFlash)
    (; numphases, max_steps, population_size, time_limit, seed, stagnation_evals, stagnation_tol, backend, verbose, logspace, equilibrium) = method
    model = __tpflash_cache_model(model,p,T,n,equilibrium)
    numspecies = length(model)
    TT = Base.promote_typeof(p, T, first(n))
    x = zeros(TT, numphases, numspecies)
    nvals = zeros(TT, numphases, numspecies)
    volumes = zeros(TT, numphases)

    # Minimize Gibbs energy
    dim = numspecies * (numphases - 1)
    lb_scalar, ub_scalar = if logspace
        (log(4eps(TT)), zero(TT))
    else
        (zero(TT), one(TT))
    end
    lb = fill(Float64(lb_scalar), dim)
    ub = fill(Float64(ub_scalar), dim)

    algo = if backend === :sass
        Solvers.SASS(population_size, max_steps, lb, ub; seed, stagnation_evals, stagnation_tol)
    else
        throw(DomainError(backend, "unknown DETPFlash backend (expected :sass)"))
    end
    deadline_ns = isfinite(time_limit) ? time_ns() + floor(Int, 1e9time_limit) : typemax(Int)
    while !Solvers.isdone(algo)
        dividers_flat = Solvers.ask!(algo)
        y = Obj_de_tp_flash(model, p, T, n, dividers_flat, numphases, x, nvals, volumes, logspace, equilibrium)
        Solvers.tell!(algo, y)
        time_ns() >= deadline_ns && break
    end
    best_u, g = Solvers.best(algo)
    # Refresh cache for the best solution (the last evaluated point might not be the best).
    Obj_de_tp_flash(model, p, T, n, copy(best_u), numphases, x, nvals, volumes, logspace, equilibrium)

    #Initialize arrays xij and nvalsij,
    #where i in 1..numphases, j in 1..numspecies
    #xij is mole fraction of j in phase i.
    #nvals is mole numbers of j in phase i.
    dividers = reshape(best_u, (numphases - 1, numspecies))
    if logspace
        dividers .= exp.(dividers)
    end
    partition!(dividers,n,x,nvals)

    comps = [vec(x[i,:]) for i in 1:numphases]
    βi = [sum(@view(nvals[i,:])) for i in 1:numphases]
    for i in 1:numphases
        if iszero(volumes[i]) && model isa PTFlashWrapper
            #we suppose liquid volume, evaluate here
            volumes[i] = volume(model,p,T,comps[i],phase = :l)
        end
    end
    return FlashResult(comps, βi, volumes, FlashData(p,T,g))
end
"""
    Obj_de_tp_flash(model,p,T,z,dividers,numphases,vcache,logspace = false)

Function to calculate Gibbs energy for given partition of moles between phases.
This is a little tricky.

We must find a way of uniquely converting a vector of numbers,
each in (0, 1), to a partition. We must be careful that
the mapping is 1-to-1, not many-to-1, as if many inputs
map to the same physical state in a redundant way, there
will be multiple global optima, and the global optimization
will perform poorly.

Our approach is to specify (numphases-1) numbers in (0,1) for
each species. We then scale these numbers systematically in order to partition
the species between the phases. Each set of (numphases - 1) numbers
will result in a unique partition of the species into the numphases
phases.
vcache stores the current volumes for each phase
"""
function Obj_de_tp_flash(model,p,T,n,dividers,numphases,x,nvals,vcache,logspace = false,equilibrium = :auto)
    # NOTE: avoid `Base.promote_typeof(p, T, first(n))` here; this function is slow
    # `x/nvals/volumes` are allocated in `tp_flash_impl` using the promoted type already, so we can reuse the cache eltype.
    TT = eltype(nvals)
    _0 = zero(TT)
    numspecies = length(n)
    bignum = TT(1e300)
    if logspace
        dividers .= exp.(dividers)
    end
    dividers = reshape(dividers, (numphases - 1, numspecies))
    #Initialize arrays xij and nvalsij,
    #where i in 1..numphases, j in 1..numspecies
    #xij is mole fraction of j in phase i.
    #nvals is mole numbers of j in phase i.
    #x = zeros(TT,numphases, numspecies)
    #nvals = zeros(TT,numphases, numspecies)
    #Calculate partition of species into phases
    partition!(dividers,n,x,nvals)
    #Calculate Overall Gibbs energy (J)
    #If any errors are encountered, return a big number, ensuring point is discarded
    #by DE Algorithm
    G = _0
    for i ∈ 1:numphases
        ni = @view(nvals[i, :])
        gi,vi = __eval_G_DETPFlash(model,p,T,ni,equilibrium)
        vcache[i] = vi
        G += gi
        #calling with PTn calls the internal volume solver
        #if it returns an error, is a bug in our part.
    end
    if logspace
        dividers .= log.(dividers)
    end
    R̄ = Rgas(model)
    return ifelse(isnan(G),bignum,G/R̄/T)
end

#indirection to allow overloading this evaluation in activity models
function __eval_G_DETPFlash(model::EoSModel,p,T,ni,equilibrium)
    phase = is_lle(equilibrium) ? :liquid : :unknown
    vi = volume(model,p,T,ni;phase = phase)
    g = VT_gibbs_free_energy(model, vi, T, ni)
    return g,vi
end

numphases(method::DETPFlash) = method.numphases

export DETPFlash
