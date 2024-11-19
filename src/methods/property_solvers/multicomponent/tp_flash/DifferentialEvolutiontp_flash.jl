#=
Original code by Thomas Moore
@denbigh
included in https://github.com/ClapeyronThermo/Clapeyron.jl/pull/56
=#
"""
    DETPFlash(; numphases = 2,
    max_steps = 1e4*(numphases-1),
    population_size =20,
    time_limit = Inf,
    verbose = false,
    logspace = false,
    equilibrium = :auto)

Method to solve non-reactive multicomponent flash problem by finding global minimum of Gibbs Free Energy via Differential Evolution.

User must assume a number of phases, `numphases`. If true number of phases is smaller than numphases, model should predict either (a) identical composition in two or more phases, or (b) one phase with negligible total number of moles. If true number of phases is larger than numphases, a thermodynamically unstable solution will be predicted.

The optimizer will stop at `max_steps` evaluations or at `time_limit` seconds

The `equilibrium` keyword allows to restrict the search of phases to just liquid-liquid equilibria (`equilibrium = :lle`). the default searches for liquid and gas phases.
"""
Base.@kwdef struct DETPFlash <: TPFlashMethod
    numphases::Int = 2
    max_steps::Int = 1e4*(numphases-1)
    population_size::Int = 50
    time_limit::Float64 = Inf
    verbose::Bool = false
    logspace::Bool = false
    equilibrium::Symbol = :auto
end

index_reduction(flash::DETPFlash,z) = flash

#z is the original feed composition, x is a matrix with molar fractions, n is a matrix with molar amounts
function partition!(dividers,n,x,nvals)
    numphases, numspecies = size(x)
    @inbounds for j = 1:numspecies
        nvals[1,j] = dividers[1,j]*n[j]
        for i = 2:numphases - 1
            Σ = sum(@view(nvals[1:(i-1), j]))
            nvals[i,j] = dividers[i,j]*(n[j] - Σ)
        end
        Σphases = sum(@view(nvals[1:(numphases-1), j]))
        nvals[numphases, j] = n[j] - Σphases
    end
#Calculate mole fractions xij
for i = 1:numphases
    ni = @view(nvals[i, :])
    invn = 1/sum(ni)
    xi = @view(x[i, :])
    xi  .= ni
    xi .*= invn
end
end

function tp_flash_impl(model::EoSModel, p, T, n, method::DETPFlash)
    model = __tpflash_cache_model(model,p,T,n,method.equilibrium)
    numspecies = length(model)
    TYPE = typeof(p+T+first(n))
    numphases = method.numphases
    x = zeros(TYPE,numphases, numspecies)
    nvals = zeros(TYPE,numphases, numspecies)
    logspace = method.logspace
    volumes = zeros(TYPE,numphases)
    GibbsFreeEnergy(dividers) = Obj_de_tp_flash(model,p,T,n,dividers,numphases,x,nvals,volumes,logspace,method.equilibrium)
    #Minimize Gibbs Free Energy

    #=
    options = Metaheuristics.Options(time_limit = method.time_limit,iterations = method.max_steps,seed = UInt(373))
    algorithm = Metaheuristics.WOA(N=method.population_size,options=options)
    bounds = vcat(zeros(TYPE,(1,numspecies*(numphases-1))),ones(TYPE,1,numspecies*(numphases-1)))
    result = Metaheuristics.optimize(GibbsFreeEnergy,bounds,algorithm)

    dividers = reshape(Metaheuristics.minimizer(result),
            (numphases - 1, numspecies))
    best_f = Metaheuristics.minimum(result)

    =#
    if logspace
        bounds = (log(4*eps(TYPE)),zero(TYPE))
    else
        bounds = (zero(TYPE),one(TYPE))
    end
    result = BlackBoxOptim.bboptimize(GibbsFreeEnergy;
        SearchRange = bounds,
        NumDimensions = numspecies*(numphases-1),
        MaxSteps=method.max_steps,
        PopulationSize = method.population_size,
        MaxTime = method.time_limit,
        TraceMode = ifelse(method.verbose,:verbose,:silent))

    dividers = reshape(BlackBoxOptim.best_candidate(result),
            (numphases - 1, numspecies))
    g = BlackBoxOptim.best_fitness(result)
    #Initialize arrays xij and nvalsij,
    #where i in 1..numphases, j in 1..numspecies
    #xij is mole fraction of j in phase i.
    #nvals is mole numbers of j in phase i.
    if logspace
        dividers .= exp.(dividers)
    end
    partition!(dividers,n,x,nvals)

    comps = [vec(x[i,:]) for i in 1:numphases]
    βi = [sum(@view(nvals[i,:])) for i in 1:numphases]
    return comps, βi, volumes, PTFlashData(p,T,g)
end
"""
    Obj_de_tp_flash(model,p,T,z,dividers,numphases,vcache,logspace = false)

Function to calculate Gibbs Free Energy for given partition of moles between phases.
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
    _0 = zero(p+T+first(n))
    numspecies = length(n)
    TYPE  = typeof(_0)
    bignum = TYPE(1e300)
    if logspace
        dividers .= exp.(dividers)
    end
    dividers = reshape(dividers,
        (numphases - 1, numspecies))
    #Initialize arrays xij and nvalsij,
    #where i in 1..numphases, j in 1..numspecies
    #xij is mole fraction of j in phase i.
    #nvals is mole numbers of j in phase i.
    #x = zeros(TYPE,numphases, numspecies)
    #nvals = zeros(TYPE,numphases, numspecies)
    #Calculate partition of species into phases
    partition!(dividers,n,x,nvals)
    #Calculate Overall Gibbs Free Energy (J)
    #If any errors are encountered, return a big number, ensuring point is discarded
    #by DE Algorithm
    G = _0
    for i ∈ 1:numphases
            xi = @view(nvals[i, :])
            gi,vi = __eval_G_DETPFlash(model,p,T,xi,equilibrium)
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
function __eval_G_DETPFlash(model::EoSModel,p,T,xi,equilibrium)
    phase = is_lle(equilibrium) ? :liquid : :unknown
    vi = volume(model,p,T,xi;phase = phase)
    g = VT_gibbs_free_energy(model, vi, T, xi)
    return g,vi
end

numphases(method::DETPFlash) = method.numphases

export DETPFlash