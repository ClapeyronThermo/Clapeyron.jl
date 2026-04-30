using Dates
using Serialization
using Printf

include(joinpath(@__DIR__, "tp_flash_benchmark_cases.jl"))
include(joinpath(@__DIR__, "tp_flash_benchmark.jl"))
include(joinpath(@__DIR__, "tp_flash_benchmark_backend.jl"))

function _sanitize_filename(s::String)
    return replace(s, r"[^A-Za-z0-9._-]+" => "_")
end

function _format_failure_tags(tags::Vector{Symbol})
    isempty(tags) && return "UNKNOWN_FAIL"
    return join(String.(tags), ", ")
end

function _git_head()
    try
        return readchomp(`git rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function tpflash_benchmark_stage1(;
    algo::AbstractString,
    backend::Symbol=:rdex,
    outdir::AbstractString=joinpath(@__DIR__, "results"),
    seeds=collect(1:10),
    population_size::Int=50,
    fe_factor::Int=10_000,
    time_limit::Float64=Inf,
    stagnation_evals::Int=0,
    stagnation_tol::Float64=0.0,
    atol_tgt::Float64=1e-14,
    rtol_tgt::Float64=1e-14,
    tau_x::Float64=1e-8,
    tau_beta::Float64=1e-8,
    beta_min::Float64=1e-12,
    max_e_x::Float64=0.01,
    max_e_beta::Float64=0.01,
    filter::Union{Nothing,Regex,AbstractString}=nothing,
    save_path::Union{Nothing,AbstractString}=nothing,
    verbose::Bool=true,
    show_fail_details::Bool=false,
)
    mkpath(outdir)

    seed_list = collect(Int, seeds)
    cfg = TPFlashBenchmarkConfig(;
        backend=backend,
        seeds=seed_list,
        population_size=population_size,
        fe_factor=fe_factor,
        time_limit=time_limit,
        stagnation_evals=stagnation_evals,
        stagnation_tol=stagnation_tol,
        atol_tgt=atol_tgt,
        rtol_tgt=rtol_tgt,
        tau_x=tau_x,
        tau_beta=tau_beta,
        beta_min=beta_min,
        max_e_x=max_e_x,
        max_e_beta=max_e_beta,
    )

    cases = tpflash_benchmark_cases()
    if filter !== nothing
        rx = filter isa Regex ? filter : Regex(String(filter))
        cases = [c for c in cases if occursin(rx, c.id)]
    end

    if verbose
        planned_save_path = save_path === nothing ? joinpath(outdir, _sanitize_filename(String(algo)) * ".ser") : String(save_path)
        println("TPFlash benchmark stage 1")
        println("- algo: ", algo)
        println("- backend: ", backend)
        println("- git:  ", _git_head())
        println("- out:   ", outdir)
        println("- save:  ", planned_save_path)
        println("- cases: ", length(cases))
        println("- seeds: ", length(seed_list), " ", seed_list)
        println("- population_size: ", population_size)
        println("- fe_factor:       ", fe_factor)
        println("- time_limit:      ", time_limit)
        println("- stagnation_evals:", stagnation_evals)
        println("- stagnation_tol:  ", stagnation_tol)
        println("- atol_tgt:        ", atol_tgt)
        println("- rtol_tgt:        ", rtol_tgt)
        println("- tau_x:           ", tau_x)
        println("- tau_beta:        ", tau_beta)
        println("- beta_min:        ", beta_min)
        println("- max_e_x:         ", max_e_x)
        println("- max_e_beta:      ", max_e_beta)
        println("- filter:          ", filter === nothing ? "nothing" : String(filter))
        println("- verbose:         ", verbose)
        println("- show_fail_details: ", show_fail_details)
    end

    results = TPFlashAlgoResults(String(algo), now(), cfg, Dict{String,Vector{TPFlashRun}}())

    for (ci, case) in enumerate(cases)
        verbose && println("\n[$ci/$(length(cases))] case ", case.id)
        ref = compute_reference(case, cfg)

        runs = TPFlashRun[]
        for (si, seed) in enumerate(seed_list)
            verbose && print("  run $(si)/$(length(seed_list)) seed=$seed ... ")
            run = run_case_with_optimizer(case, ref, seed, cfg)
            push!(runs, run)
            if verbose
                if run.success
                    println("success")
                elseif show_fail_details
                    @printf("%s Δ=%0.3e (g=%0.16g tgt=%0.16g) e_x=%0.3e e_beta=%0.3e nfes=%d\n",
                        _format_failure_tags(run.failure_tags),
                        run.Δ_final,
                        run.g_best,
                        run.tgt,
                        run.e_x,
                        run.e_beta,
                        run.nfes,
                    )
                else
                    println(_format_failure_tags(run.failure_tags))
                end
            end
        end

        results.case_runs[case.id] = runs
    end

    save_path = save_path === nothing ? joinpath(outdir, _sanitize_filename(String(algo)) * ".ser") : String(save_path)
    serialize(save_path, results)
    if verbose
        println()
        println("Wrote: ", save_path)
    end

    return results, save_path
end

function _argval(key::String; default=nothing)
    prefix = "--$key="
    for a in ARGS
        startswith(a, prefix) && return a[length(prefix)+1:end]
    end
    return default
end

function _parse_int(key; default)
    v = _argval(String(key); default=nothing)
    v === nothing && return default
    return parse(Int, v)
end

function _parse_float(key; default)
    v = _argval(String(key); default=nothing)
    v === nothing && return default
    v == "Inf" && return Inf
    return parse(Float64, v)
end

function _parse_seeds()
    v = _argval("seeds"; default=nothing)
    if v === nothing
        n = _parse_int("n"; default=10)
        start = _parse_int("seed-start"; default=1)
        return collect(start:(start+n-1))
    end
    return [parse(Int, s) for s in split(v, ",") if !isempty(s)]
end

function main()
    algo = _argval("algo"; default=nothing)
    algo === nothing && error("missing required arg: --algo=<name>")
    backend = Symbol(_argval("backend"; default="rdex"))

    outdir = _argval("outdir"; default=joinpath(@__DIR__, "results"))
    seeds = _parse_seeds()
    pop = _parse_int("pop"; default=50)
    fe_factor = _parse_int("fe-factor"; default=10_000)
    time_limit = _parse_float("time-limit"; default=Inf)
    stagnation_evals = _parse_int("stagnation-evals"; default=0)
    stagnation_tol = _parse_float("stagnation-tol"; default=0.0)
    atol_tgt = _parse_float("atol-tgt"; default=1e-14)
    rtol_tgt = _parse_float("rtol-tgt"; default=1e-14)
    tau_x = _parse_float("tau-x"; default=1e-8)
    tau_beta = _parse_float("tau-beta"; default=1e-8)
    beta_min = _parse_float("beta-min"; default=1e-12)
    max_e_x = _parse_float("max-e-x"; default=0.01)
    max_e_beta = _parse_float("max-e-beta"; default=0.01)
    filter = _argval("filter"; default=nothing)
    tpflash_benchmark_stage1(;
        algo=algo,
        backend=backend,
        outdir=outdir,
        seeds=seeds,
        population_size=pop,
        fe_factor=fe_factor,
        time_limit=time_limit,
        stagnation_evals=stagnation_evals,
        stagnation_tol=stagnation_tol,
        atol_tgt=atol_tgt,
        rtol_tgt=rtol_tgt,
        tau_x=tau_x,
        tau_beta=tau_beta,
        beta_min=beta_min,
        max_e_x=max_e_x,
        max_e_beta=max_e_beta,
        filter=filter,
        verbose=true,
    )
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
