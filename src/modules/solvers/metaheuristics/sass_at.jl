mutable struct Arch
    max_npop::Int
    x::Matrix{Float64}
    y::Vector{Float64}
end

struct SASSConfig
    max_npop   :: Int
    min_npop   :: Int
    max_nfes   :: Int
    stagnation_evals::Int
    stagnation_tol::Float64
    pbestratio :: Float64
    Ar         :: Float64
    Ms         :: Int
    dim        :: Int
    lb         :: Vector{Float64}
    ub         :: Vector{Float64}
    seed       :: Int
end

mutable struct SASSState
    gen      :: Int
    npop     :: Int
    nfes     :: Int
    stagnation_count::Int

    rng        :: MersenneTwister

    # populations
    x        :: Matrix{Float64}
    x_new    :: Matrix{Float64}
    x_shadow :: Matrix{Float64}
    y        :: Vector{Float64}         # previous generation fitness (baseline)
    curr_y2  :: Vector{Float64}         # current generation collected fitness

    # memory and archive
    arch   :: Arch
    ci     :: Vector{Float64}
    rd     :: Vector{Float64}
    mem_c̄ :: Vector{Float64}
    mem_rd :: Vector{Float64}

    # best tracking
    best_x::Vector{Float64}
    best_y_current::Float64
    best_y_history::Vector{Float64}

    # ask–tell bookkeeping
    is_pending  :: Bool
    pending_col :: Int
    pending_x   :: Vector{Float64}
    eval_col    :: Int                    # next column to evaluate (1..npop)
end

mutable struct SASSWorkspace
    perm_buf::Vector{Int}
end

mutable struct SASS
    config  :: SASSConfig
    state   :: SASSState
    buffers :: SASSWorkspace
end

function SASS(
    max_npop::Integer,
    max_nfes::Integer,
    lb::AbstractVector,
    ub::AbstractVector;
    seed::Integer = 1,
    stagnation_evals::Integer = 0,
    stagnation_tol::Real = 0.0,
)
    min_npop = 4
    pbestratio = 0.11
    rd0 = 0.9
    c̄0 = 0.7
    Ar = 1.4
    Ms = 5
    dim = length(lb)

    max_npop = Int(max_npop)
    max_nfes = Int(max_nfes)
    stagnation_evals = Int(stagnation_evals)
    stagnation_tol = Float64(stagnation_tol)
    max_npop >= 4 || throw(DomainError(max_npop, "Population size must be at least 4"))
    max_nfes >= max_npop || throw(DomainError(max_nfes, "max_nfes must be >= population size"))
    stagnation_evals >= 0 || throw(DomainError(stagnation_evals, "stagnation_evals must be >= 0 (0 disables early stopping)"))
    stagnation_tol >= 0.0 || throw(DomainError(stagnation_tol, "stagnation_tol must be >= 0"))

    cfg = SASSConfig(max_npop, min_npop, max_nfes,
        stagnation_evals, stagnation_tol,
        pbestratio, Ar, Ms,
        dim, collect(Float64, lb), collect(Float64, ub), Int(seed))

    npop = max_npop
    gen = 1
    nfes = 0
    stagnation_count = 0

    rng = MersenneTwister(cfg.seed)

    x = get_population((dim, npop), cfg.lb, cfg.ub, rng)
    x_new = copy(x)
    x_shadow = get_population((dim, npop), cfg.lb, cfg.ub, rng)
    y = Vector{Float64}(undef, 0)                    # baseline filled after gen 1
    curr_y2 = fill(NaN, npop)
    best_x = Vector{Float64}(undef, dim)
    best_y_current = Inf
    best_y_history = sizehint!(Float64[], max_nfes ÷ max(1, max_npop) * 5)

    arch = Arch(round(Int, Ar * npop), zeros(dim, 0), Float64[])
    ci = Float64[]
    rd = Float64[]
    mem_rd = fill(rd0, Ms)
    mem_c̄ = fill(c̄0, Ms)

    state = SASSState(gen, npop, nfes, stagnation_count,
        rng,
        x, x_new, x_shadow, y, curr_y2,
        arch, ci, rd, mem_c̄, mem_rd,
        best_x, best_y_current, best_y_history,
        false, 1, Vector{Float64}(undef, dim), 1)

    buffers = SASSWorkspace(zeros(Int, npop))
    return SASS(cfg, state, buffers)
end

@inline function _gen_r1r2r3(rng::AbstractRNG, npop::Integer, popAllSize::Integer)
    r0 = 1:npop
    r3 = randperm(rng, npop)
    i = 1
    while any(r3 .== r0)
        randperm!(rng, r3)
        i += 1
        i > 1000 && error("Cannot generate r3 in 1000 iterations")
    end
    r1 = Vector{Int}(undef, npop)
    for i = 1:npop
        r1[i] = rand(rng, 1:npop)
        retries = 0
        @inbounds while r1[i] == r0[i] || r1[i] == r3[i]
            r1[i] = rand(rng, 1:npop)
            retries += 1
            retries > 1000 && error("Cannot generate r1[$i] in 1000 iterations")
        end
    end
    r2 = Vector{Int}(undef, npop)
    for i = 1:npop
        r2[i] = rand(rng, 1:popAllSize)
        retries = 0
        @inbounds while r2[i] == r0[i] || r2[i] == r1[i] || r2[i] == r3[i]
            r2[i] = rand(rng, 1:popAllSize)
            retries += 1
            retries > 1000 && error("Cannot generate r2[$i] in 1000 retries")
        end
    end
    return r1, r2, r3
end

"Orthogonal matrix approx"
function rand_orth_mat(rng::AbstractRNG, n::Integer, t::Real)
    R = Matrix{Float64}(I, n, n)
    l = randperm(rng, n)
    for i = 1:2:n-1
        R[l[i], l[i]] = sin(t)
        R[l[i+1], l[i+1]] = sin(t)
        R[l[i], l[i+1]] = cos(t)
        R[l[i+1], l[i]] = -cos(t)
    end
    return R
end

function upd_arch!(rng::AbstractRNG, arch::Arch, x::AbstractMatrix, funvalue::AbstractVector)
    arch.max_npop == 0 && return
    size(x, 2) == length(funvalue) || throw(DimensionMismatch("size(x, 2) != length(funvalue)"))
    pop_all = [arch.x x]
    pop_all = unique(pop_all, dims=2)
    actualPop = size(pop_all, 2)
    if actualPop <= arch.max_npop
        arch.x = pop_all
    else
        rndPos = randperm(rng, actualPop)[1:arch.max_npop]
        arch.x = pop_all[:, rndPos]
    end
    nothing
end

@inline function _emit_next_candidate!(sa::SASS)
    s = sa.state
    if s.is_pending
        return s.pending_x
    end
    # find next unevaluated index in curr_y2
    i = findfirst(isnan, s.curr_y2)
    if i === nothing
        return nothing
    end
    s.pending_col = i
    s.pending_x .= @view s.x_new[:, i]
    s.is_pending = true
    return s.pending_x
end

function _finalize_generation!(sa::SASS)
    c = sa.config
    s = sa.state
    b = sa.buffers
    rng = s.rng

    i = s.gen
    npop = s.npop
    dim = c.dim
    lb, ub = c.lb, c.ub
    Ms = c.Ms
    Ar = c.Ar
    pbestratio = c.pbestratio

    y₂ = s.curr_y2
    push!(s.best_y_history, s.best_y_current)
    # `nfes` is incremented in `tell!` (one evaluation per tell); do not accumulate here.

    if i == 1
        s.y = copy(y₂)
    else
        gdidx = y₂ .< s.y
        gdnum = count(gdidx)
        if gdnum > 0
            GdRD = s.rd[gdidx]
            GdCI = s.ci[gdidx]
            diffVal = s.y[gdidx] .- y₂[gdidx]
            upd_arch!(rng, s.arch, s.x[:, gdidx], s.y[gdidx])
            for j in findall(gdidx)
                s.x_shadow[:, j] .= @view s.x[:, j]
                s.x[:, j] .= @view s.x_new[:, j]
                s.y[j] = y₂[j]
            end
            mem_I = ((i - 2) % Ms) + 1
            diffVal ./= sum(diffVal)
            s.mem_c̄[mem_I] = (diffVal ⋅ (GdCI .^ 2)) / (diffVal ⋅ GdCI)
            if maximum(GdRD) == 0 || s.mem_rd[mem_I] == -1
                s.mem_rd[mem_I] = -1
            else
                s.mem_rd[mem_I] = (diffVal ⋅ (GdRD .^ 2)) / (diffVal ⋅ GdRD)
            end
        end

        plan_popSize = round(Int, c.min_npop * (s.nfes / c.max_nfes) + c.max_npop * (1 - s.nfes / c.max_nfes))
        plan_popSize = max(plan_popSize, c.min_npop)
        if npop > plan_popSize
            sorted_indices = sortperm!(b.perm_buf, s.y)
            to_remove = sort!(view(sorted_indices, plan_popSize+1:npop))
            deleteat!(s.y, to_remove)
            resize!(b.perm_buf, plan_popSize)
            survivors = sort!(b.perm_buf)
            s.x = s.x[:, survivors]
            s.x_shadow = s.x_shadow[:, survivors]
            npop = plan_popSize
            s.npop = npop
            s.arch.max_npop = round(Int, Ar * npop)
            actualPop = size(s.arch.x, 2)
            if actualPop > s.arch.max_npop
                rndPos = randperm(rng, actualPop)[1:s.arch.max_npop]
                s.arch.x = s.arch.x[:, rndPos]
            end
        end
    end

    Base.permutecols!!(s.x_shadow, randperm!(rng, sa.buffers.perm_buf))
    mem_rand_index = rand!(rng, sa.buffers.perm_buf, 1:Ms)
    MUci = s.mem_c̄[mem_rand_index]
    MUrd = s.mem_rd[mem_rand_index]
    s.rd = (0.1 * randn(rng)) .+ MUrd
    s.rd[MUrd.==-1] .= 0
    clamp!(s.rd, 0, 1)
    r_ci = rand(rng, s.npop)
    s.ci = @. MUci + 0.1 * tan(π * (r_ci - 0.5))
    Pos = s.ci .<= 0
    while any(Pos)
        npos = count(identity, Pos)
        s.ci[Pos] .= MUci[Pos] .+ 0.1tan.(π .* (rand(rng, npos) .- 0.5))
        Pos .= s.ci .<= 0
    end
    clamp!(s.ci, 0, 1)

    pop_all = [s.x s.arch.x]
    (r1, r2, r3) = _gen_r1r2r3(rng, s.npop, size(pop_all, 2))
    pNP = max(round(Int, pbestratio * s.npop), 2)
    sorted_indices = sortperm!(sa.buffers.perm_buf, s.y)
    elite_idx = rand(rng, view(sorted_indices, 1:pNP), s.npop)
    pbest = s.x[:, elite_idx]
    ks = sorted_indices[1:floor(Int, s.npop / 2)]
    pbest[:, ks] .= @view s.x[:, r3[ks]]
    popA = s.x_shadow[:, r1]
    popA[:, ks] .= @view pop_all[:, r2[ks]]
    z = @. pbest + s.x[:, r1] - popA - s.x
    A = rand_orth_mat(rng, dim, 1e-12)
    mul!(pbest, A', z)
    copy!(z, pbest)
    Ur = fill!(pbest, 0.0)
    for i = 1:s.npop
        r = rand(rng, 1:dim)
        idx = (r - 1) * s.npop + i
        Ur[idx] = z[idx]
    end
    for j = 1:s.npop, i = 1:dim
        if rand(rng) < s.rd[j]
            Ur[i, j] = z[i, j]
        end
    end
    z .= s.ci' .* Ur
    step = pbest
    mul!(step, A, z)
    confine!(step, s.x, lb, ub, 0.5)
    s.x_new = z .= s.x .+ step

    s.gen += 1
    s.curr_y2 = fill(NaN, s.npop)
    s.eval_col = 1
    s.is_pending = false
    return nothing
end

function ask!(sa::SASS)
    s = sa.state
    if s.is_pending
        return s.pending_x
    end
    # try emit pending within current generation
    x = _emit_next_candidate!(sa)
    if x !== nothing
        return x
    end
    # current generation finished → build next generation and emit first
    _finalize_generation!(sa)
    return _emit_next_candidate!(sa)
end

function tell!(sa::SASS, y::Real)
    s = sa.state
    c = sa.config
    s.is_pending || error("no pending candidate to tell!")
    y = Float64(y)
    i = s.pending_col
    s.curr_y2[i] = y

    old_best_y = s.best_y_current
    if y < s.best_y_current
        s.best_y_current = y
        s.best_x .= s.pending_x
    end
    if c.stagnation_evals > 0
        if (y + c.stagnation_tol) < old_best_y
            s.stagnation_count = 0
        else
            s.stagnation_count += 1
        end
    end

    s.nfes += 1
    s.is_pending = false
    # if still have unevaluated individuals, continue; otherwise generation will be finalized on next ask!
    return nothing
end

function isdone(sa::SASS)
    s = sa.state
    c = sa.config
    stagnating = (c.stagnation_evals > 0) && (s.stagnation_count >= c.stagnation_evals)
    return (s.nfes >= c.max_nfes) || stagnating
end

best(sa::SASS) = (sa.state.best_x, sa.state.best_y_current)
