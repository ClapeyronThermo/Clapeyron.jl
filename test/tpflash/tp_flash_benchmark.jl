using Clapeyron
using Serialization
using Dates

"""
Benchmark configuration (defaults follow `tp_flash_definitions_and_eval_criteria.tex`).
"""
Base.@kwdef struct TPFlashBenchmarkConfig
    backend::Symbol = :rdex              # :rdex (ask–tell) or :bbo (BlackBoxOptim)
    seeds::Vector{Int} = collect(1:10)
    fe_factor::Int = 10_000                 # FE_max = fe_factor * D
    population_size::Int = 50
    time_limit::Float64 = Inf
    stagnation_evals::Int = 0               # passed to optimizer backend; 0 disables
    stagnation_tol::Float64 = 0.0

    atol_tgt::Float64 = 1e-14
    rtol_tgt::Float64 = 1e-14

    tau_x::Float64 = 1e-8
    tau_beta::Float64 = 1e-8
    beta_min::Float64 = 1e-12

    max_e_x::Float64 = 0.01
    max_e_beta::Float64 = 0.01
end

struct TPFlashReference
    g_star::Float64
    tgt::Float64
    x_star::Matrix{Float64}   # (numphases, numspecies)
    beta_star::Vector{Float64}# length numphases, normalized
end

struct TPFlashRun
    seed::Int
    fe_max::Int
    nfes::Int
    fe_star::Union{Int,Nothing}
    g_best::Float64
    tgt::Float64
    ev_final::Float64
    success::Bool
    failure_tags::Vector{Symbol}  # may include multiple tags (e.g. TARGET_FAIL, X_FAIL, BETA_FAIL)
    e_x::Float64
    e_beta::Float64
    x_best::Matrix{Float64}
    beta_best::Vector{Float64}
    u_best::Vector{Float64}
end

struct TPFlashAlgoResults
    algo::String
    created_at::DateTime
    config::TPFlashBenchmarkConfig
    case_runs::Dict{String,Vector{TPFlashRun}}
end

function _tgt_from_gstar(gstar::Float64; atol::Float64, rtol::Float64)
    tol = max(100 * eps(gstar), rtol * abs(gstar) + atol)
    return gstar + tol
end

function _normalize_beta(beta::AbstractVector{<:Real})
    s = sum(beta)
    if !(s > 0)
        return fill(NaN, length(beta))
    end
    return collect(Float64, beta) ./ Float64(s)
end

const _PERMS_1 = ((1,),)
const _PERMS_2 = ((1, 2), (2, 1))
const _PERMS_3 = (
    (1, 2, 3),
    (1, 3, 2),
    (2, 1, 3),
    (2, 3, 1),
    (3, 1, 2),
    (3, 2, 1),
)

function _perms(k::Int)
    if k == 1
        return _PERMS_1
    elseif k == 2
        return _PERMS_2
    elseif k == 3
        return _PERMS_3
    else
        error("phase matching only implemented for k<=3 (got k=$k)")
    end
end

function _phase_match_errors(
    x::AbstractMatrix{<:Real},
    beta::AbstractVector{<:Real},
    xref::AbstractMatrix{<:Real},
    betaref::AbstractVector{<:Real};
    beta_min::Float64,
)
    k = size(x, 1)
    size(xref, 1) == k || error("phase count mismatch (x has $k, xref has $(size(xref, 1)))")
    length(beta) == k || error("beta length mismatch")
    length(betaref) == k || error("beta_ref length mismatch")

    best_perm = collect(1:k)
    best_e_beta = Inf
    best_e_x = Inf

    for perm in _perms(k)
        e_beta = 0.0
        e_x = 0.0
        @inbounds for i = 1:k
            bi = Float64(beta[i])
            bj = Float64(betaref[perm[i]])
            active = max(bi, bj) > beta_min
            if active
                e_beta = max(e_beta, abs(bi - bj))
                e_x = max(e_x, Clapeyron.dnorm(view(x, i, :), view(xref, perm[i], :), 2))
            end
        end

        if (e_beta < best_e_beta) || (e_beta == best_e_beta && e_x < best_e_x)
            best_perm = collect(perm)
            best_e_beta = e_beta
            best_e_x = e_x
        end
    end

    return best_perm, best_e_x, best_e_beta
end

function _validity_gate(
    g::Float64,
    x::AbstractMatrix{<:Real},
    beta::AbstractVector{<:Real};
    tau_x::Float64,
    tau_beta::Float64,
    beta_min::Float64,
)
    if !isfinite(g)
        return false, :NUMERIC_FAIL
    end
    any(!isfinite, x) && return false, :NUMERIC_FAIL
    any(!isfinite, beta) && return false, :NUMERIC_FAIL

    sβ = sum(beta)
    if any(b -> b < -tau_beta, beta) || !(abs(sβ - 1.0) <= tau_beta)
        return false, :CONSTRAINT_FAIL
    end

    k, _ = size(x)
    @inbounds for i = 1:k
        if beta[i] > beta_min
            row = view(x, i, :)
            if any(v -> v < -tau_x, row)
                return false, :CONSTRAINT_FAIL
            end
            sx = sum(row)
            if !(abs(sx - 1.0) <= tau_x)
                return false, :CONSTRAINT_FAIL
            end
        end
    end

    return true, :OK
end

"""
Return the hardcoded reference (g*, x*, β*) for one benchmark case.

Reference values are stored as numeric literals in `test/tpflash/tp_flash_benchmark_cases.jl`
so the benchmark does not depend on online reference solving.
"""
function compute_reference(case, cfg::TPFlashBenchmarkConfig)
    (case.g⁺ === nothing || case.x⁺ === nothing || case.β⁺ === nothing) && error("missing reference literals for case $(case.id)")
    g = Float64(case.g⁺)
    x = Matrix{Float64}(case.x⁺)
    beta = Vector{Float64}(case.β⁺)

    if size(x, 1) != case.numphases
        x_pad = zeros(Float64, case.numphases, size(x, 2))
        beta_pad = zeros(Float64, case.numphases)
        ncopy = min(case.numphases, size(x, 1))
        x_pad[1:ncopy, :] .= x[1:ncopy, :]
        beta_pad[1:ncopy] .= beta[1:ncopy]
        x = x_pad
        beta = _normalize_beta(beta_pad)
    end

    tgt = _tgt_from_gstar(g; atol=cfg.atol_tgt, rtol=cfg.rtol_tgt)
    return TPFlashReference(g, tgt, x, beta)
end

"""
Solve a TP flash benchmark case with an ask–tell optimizer over the DETPFlash objective.

Optimizer creation is delegated to `tpflash_make_optimizer` (defined in
`test/tpflash/tp_flash_benchmark_backend.jl`) so you can swap algorithms by hacking
that file without changing the benchmark harness.
"""
function run_case_with_optimizer(case, ref::TPFlashReference, seed::Int, cfg::TPFlashBenchmarkConfig)
    model0 = case.model_builder()
    model = Clapeyron.__tpflash_cache_model(model0, case.p, case.T, case.feed, case.equilibrium)

    numspecies = length(model)
    numphases = case.numphases
    dim = numspecies * (numphases - 1)
    fe_max = cfg.fe_factor * dim

    lb_scalar, ub_scalar = if case.logspace
        (log(4eps(Float64)), 0.0)
    else
        (0.0, 1.0)
    end
    lb = fill(Float64(lb_scalar), dim)
    ub = fill(Float64(ub_scalar), dim)

    x_cache = zeros(Float64, numphases, numspecies)
    nvals_cache = zeros(Float64, numphases, numspecies)
    volumes_cache = zeros(Float64, numphases)

    eval_count = Ref(0)
    fe_star_ref = Ref{Union{Int,Nothing}}(nothing)
    best_so_far = Ref(Inf)

    function objective(u)
        eval_count[] += 1
        u_eval = case.logspace ? copy(u) : u
        y = try
            Clapeyron.Obj_de_tp_flash(model, case.p, case.T, case.feed, u_eval, numphases,
                x_cache, nvals_cache, volumes_cache, case.logspace, case.equilibrium)
        catch
            Inf
        end
        y = Float64(y)
        if y < best_so_far[]
            best_so_far[] = y
        end
        if fe_star_ref[] === nothing && best_so_far[] <= ref.tgt
            fe_star_ref[] = eval_count[]
        end
        return y
    end

    best_u, best_g = tpflash_optimize(objective;
        backend=cfg.backend,
        population_size=cfg.population_size,
        fe_max=fe_max,
        lb=lb,
        ub=ub,
        seed=seed,
        time_limit=cfg.time_limit,
        stagnation_evals=cfg.stagnation_evals,
        stagnation_tol=cfg.stagnation_tol,
        verbose=false,
    )
    nfes = eval_count[]
    fe_star = fe_star_ref[]

    best_u_eval = case.logspace ? copy(best_u) : best_u
    g_best, beta, x = try
        best_g2 = Clapeyron.Obj_de_tp_flash(model, case.p, case.T, case.feed, best_u_eval, numphases, x_cache, nvals_cache, volumes_cache, case.logspace, case.equilibrium)
        g_best = Float64(best_g2)
        beta_amounts = [sum(@view(nvals_cache[i, :])) for i = 1:numphases]
        beta = _normalize_beta(beta_amounts)
        x = copy(x_cache)
        g_best, beta, x
    catch
        Inf, fill(NaN, numphases), fill(NaN, numphases, numspecies)
    end
    u_best = collect(Float64, best_u)

    valid, vcode = _validity_gate(g_best, x, beta; tau_x=cfg.tau_x, tau_beta=cfg.tau_beta, beta_min=cfg.beta_min)
    if !valid
        return TPFlashRun(seed, fe_max, nfes, fe_star, g_best, ref.tgt, g_best - ref.tgt, false, Symbol[vcode], NaN, NaN, x, beta, u_best)
    end

    _, e_x, e_beta = _phase_match_errors(x, beta, ref.x_star, ref.beta_star; beta_min=cfg.beta_min)

    tags = Symbol[]
    g_best > ref.tgt && push!(tags, :TARGET_FAIL)
    e_x > cfg.max_e_x && push!(tags, :X_FAIL)
    e_beta > cfg.max_e_beta && push!(tags, :BETA_FAIL)

    if isempty(tags)
        return TPFlashRun(seed, fe_max, nfes, fe_star, g_best, ref.tgt, g_best - ref.tgt, true, Symbol[], e_x, e_beta, x, beta, u_best)
    else
        return TPFlashRun(seed, fe_max, nfes, fe_star, g_best, ref.tgt, g_best - ref.tgt, false, tags, e_x, e_beta, x, beta, u_best)
    end
end

function _avg_rank!(ranks::Vector{Float64}, sorted_keys)
    n = length(sorted_keys)
    i = 1
    while i <= n
        j = i
        while j < n && sorted_keys[j+1] == sorted_keys[i]
            j += 1
        end
        avg = (i + j) / 2
        for k = i:j
            ranks[k] = avg
        end
        i = j + 1
    end
    return ranks
end

"""
Compute per-run ranks (1..N, where N is best) using the protocol:

- Any success outranks any failure
- Success vs success: smaller EV_final (or g_best) is better; FE* breaks ties
- Failure vs failure: smaller EV_final is better

Returns ranks aligned with `runs` order.
"""
function ranks_for_case(runs::Vector{TPFlashRun})
    n = length(runs)
    idx = collect(1:n)

    function key(r::TPFlashRun)
        if r.success
            fe = something(r.fe_star, r.nfes)
            # Sort from worst -> best (ascending key), so the best run gets the largest rank (N).
            # Successes always outrank failures. Within successes, prefer smaller EV_final (lower g_best),
            # and use FE* (time-to-target) only as a tie-breaker.
            return (1, -r.ev_final, -Float64(fe))
        else
            return (0, -r.ev_final, 0.0) # failures first; larger rank for smaller EV_final
        end
    end

    keys = map(i -> key(runs[i]), idx)
    p = sortperm(1:n; by=i -> keys[i])
    idx_sorted = idx[p]
    keys_sorted = keys[p]

    ranks_sorted = Vector{Float64}(undef, n)
    _avg_rank!(ranks_sorted, keys_sorted)

    ranks = Vector{Float64}(undef, n)
    @inbounds for pos = 1:n
        ranks[idx_sorted[pos]] = ranks_sorted[pos]
    end
    return ranks
end

function u_score_for_case(runs_by_algo::Dict{String,Vector{TPFlashRun}})
    algos = collect(keys(runs_by_algo))
    sort!(algos)
    m = length(algos)
    n = length(first(values(runs_by_algo)))
    N = m * n

    all_runs = TPFlashRun[]
    run_algo = String[]
    for a in algos
        rs = runs_by_algo[a]
        length(rs) == n || error("inconsistent seeds per algorithm for case (algo=$a)")
        append!(all_runs, rs)
        append!(run_algo, fill(a, length(rs)))
    end

    ranks = ranks_for_case(all_runs)

    sr = Dict{String,Float64}(a => 0.0 for a in algos)
    for (r, a) in zip(ranks, run_algo)
        sr[a] += r
    end

    cf = n * (n + 1) / 2
    u = Dict{String,Float64}()
    u_norm = Dict{String,Float64}()
    denom = n * (N - n)
    for a in algos
        ua = sr[a] - cf
        u[a] = ua
        u_norm[a] = denom > 0 ? ua / denom : NaN
    end

    return (; algos, ranks, run_algo, u, u_norm)
end
