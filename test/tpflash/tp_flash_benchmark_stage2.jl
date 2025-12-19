using Serialization
using Dates

include(joinpath(@__DIR__, "tp_flash_benchmark.jl"))

function _argval(key::String; default = nothing)
    prefix = "--$key="
    for a in ARGS
        startswith(a, prefix) && return a[length(prefix) + 1:end]
    end
    return default
end

function _sanitize_filename(s::String)
    return replace(s, r"[^A-Za-z0-9._-]+" => "_")
end

function _write_csv(path::AbstractString, header::Vector{String}, rows::Vector{Vector})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ","))
        for r in rows
            vals = String[]
            for v in r
                if v isa String
                    push!(vals, "\"" * replace(v, "\"" => "\"\"") * "\"")
                elseif v isa Symbol
                    push!(vals, String(v))
                else
                    push!(vals, string(v))
                end
            end
            println(io, join(vals, ","))
        end
    end
end

function tpflash_benchmark_stage2(;
    indir::AbstractString = joinpath(@__DIR__, "results"),
    outdir::AbstractString = String(indir),
    files::Union{Nothing, Vector{String}} = nothing,
    verbose::Bool = true,
)
    mkpath(outdir)

    local_files = if files === nothing
        fs = filter(f -> endswith(lowercase(f), ".ser"), readdir(indir; join = true))
        fs
    else
        [isabspath(f) ? f : joinpath(indir, f) for f in files]
    end
    isempty(local_files) && error("no .ser files found in $indir")

    alg_results = TPFlashAlgoResults[]
    for f in local_files
        obj = deserialize(f)
        if obj isa TPFlashAlgoResults
            push!(alg_results, obj)
        else
            verbose && println("Skipping non-result .ser file: ", f)
        end
    end
    isempty(alg_results) && error("no TPFlashAlgoResults found in $indir")

    algos = sort([r.algo for r in alg_results])
    algo_by_name = Dict(r.algo => r for r in alg_results)

    common_cases = nothing
    common_seeds = nothing
    for r in alg_results
        cases = sort(collect(keys(r.case_runs)))
        seeds = r.config.seeds
        if common_cases === nothing
            common_cases = cases
            common_seeds = seeds
        else
            cases == common_cases || error("case set mismatch for algo=$(r.algo)")
            seeds == common_seeds || error("seed list mismatch for algo=$(r.algo)")
        end
    end
    case_ids = common_cases::Vector{String}
    seeds = common_seeds::Vector{Int}
    n = length(seeds)

    if verbose
        println("TPFlash benchmark stage 2")
        println("- algos: ", length(algos))
        println("- cases: ", length(case_ids))
        println("- seeds: ", n)
        println("- in:    ", indir)
        println("- out:   ", outdir)
        if length(algos) == 1
            println("- note:  only one algorithm; u_norm is NaN (needs >=2 algorithms to compare)")
        end
    end

    per_case_rows = Vector{Vector}()
    summary_u_norm = Dict(a => 0.0 for a in algos)
    summary_success = Dict(a => 0 for a in algos)
    summary_total = Dict(a => 0 for a in algos)

    for case_id in case_ids
        runs_by_algo = Dict{String, Vector{TPFlashRun}}()
        for a in algos
            runs_by_algo[a] = algo_by_name[a].case_runs[case_id]
        end
        res = u_score_for_case(runs_by_algo)

        for a in algos
            u_norm = res.u_norm[a]
            summary_u_norm[a] += u_norm
            rs = runs_by_algo[a]
            summary_success[a] += count(r -> r.success, rs)
            summary_total[a] += length(rs)

            per_case_rows_push = Any[
                case_id,
                a,
                u_norm,
                res.u[a],
                count(r -> r.success, rs) / length(rs),
            ]
            push!(per_case_rows, per_case_rows_push)
        end
    end

    p = length(case_ids)
    summary_rows = Vector{Vector}()
    for a in algos
        mean_u_norm = summary_u_norm[a] / p
        sr = summary_total[a] == 0 ? 0.0 : summary_success[a] / summary_total[a]
        push!(summary_rows, Any[a, mean_u_norm, sr, summary_success[a], summary_total[a]])
    end
    sort!(summary_rows; by = r -> r[2], rev = true)

    stamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
    out_base = _sanitize_filename("tpflash_benchmark_" * stamp)
    per_case_path = joinpath(outdir, out_base * "_per_case.csv")
    summary_path = joinpath(outdir, out_base * "_summary.csv")

    _write_csv(
        per_case_path,
        ["case_id", "algo", "u_norm", "u", "success_rate"],
        per_case_rows,
    )
    _write_csv(
        summary_path,
        ["algo", "mean_u_norm", "overall_success_rate", "successes", "runs"],
        summary_rows,
    )

    if verbose
        println()
        println("Wrote:")
        println("- ", summary_path)
        println("- ", per_case_path)
    end

    return (; summary_path, per_case_path, algos, case_ids, seeds, summary_rows)
end

function main()
    indir = _argval("indir"; default = joinpath(@__DIR__, "results"))
    outdir = _argval("outdir"; default = indir)
    tpflash_benchmark_stage2(; indir = indir, outdir = outdir, verbose = true)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
