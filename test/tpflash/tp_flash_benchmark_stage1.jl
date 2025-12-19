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
                    @printf("%s ev=%0.3e (g=%0.16g tgt=%0.16g) e_x=%0.3e e_beta=%0.3e nfes=%d\n",
                        _format_failure_tags(run.failure_tags),
                        run.ev_final,
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

#=
TODO:
我刚把你 test/tpflash/results/rdex.ser 反序列化看了下，失败大致分三类（所以不是单纯 tol 的问题）：

1.差一点点的 TARGET_FAIL（基本就是 tol 太苛刻）
vlle_pcsaft_water_cyclohexane_propane_298K_1bar：g_best=-6.759674475174063，tgt=-6.759674475174976，g_best-tgt=+9.13e-13（非常接近，但仍被判 TARGET_FAIL）。
2.根本没搜到好的解（不是 tol）
vle_454_case2_rr：g_best≈580.78，tgt≈103.70（差 477）。
vle_454_case7_michelsen：g_best≈8565.53，tgt≈112.96（差 8452）。
这类需要靠算法/参数（fe_factor、population_size、以及很可能给 #454 打开 logspace=true）去改善。
X_FAIL / BETA_FAIL（目标到了但 x/β 不像参考）
3. vle_pr_*_282p2K、vle_pr_*_450K：g_best <= tgt（g_best-tgt 是负的，说明过了 target gate），但 e_beta 很大（~0.19/0.30），所以卡在 accuracy gate。它们本质上是“接近单相退化”的情形，β 往往不容易被唯一确定/优化到参考那种 0/1 分配。
为了让你下一次不需要反序列化看，我已经在 test/tpflash/tp_flash_benchmark_stage1.jl 加了开关：show_fail_details=true 时会把 ev_final(g_best-tgt)/e_x/e_beta/nfes 直接打印出来。


3. 这里的关键点是：“最小 Gibbs 能量”并不一定对应唯一的 β，尤其在你这两个 case 里其实是“单相/近单相退化”的情形。
以 DETPFlash 的目标函数 Obj_de_tp_flash 为例，它是在你假设的 numphases=2 下，把总摩尔数按 dividers 分到两相，然后算
[
G = G_1(n_1,x_1) + G_2(n_2,x_2)
]
再除以 (RT) 得到标量目标。

如果真实物理状态是单相（或两相组成几乎一样），那么最优很可能出现

两相组成几乎相同：x₁ ≈ x₂ ≈ z（你看到 x_best 对 [0.25,0.25,0.25,0.25] 的误差只有 1e-8~1e-7 量级）
此时摩尔 Gibbs 能（或化学势）在两相里几乎一样，所以把同样的混合物“分成两份”并不会改变总 (G)：
直观上就是 (G \approx n_1 \bar g(z) + n_2 \bar g(z) = (n_1+n_2)\bar g(z))，对怎么分（β）基本不敏感。
所以会出现你说的现象：g_best 已经过了 tgt（target gate 只看 g），但 β 还可以在一个很宽的范围内漂（accuracy gate 里 e_beta 就会偏大）。这里的 β⁺=[1,0] 更多是一种“约定”（把第二相视为消失），不是由目标函数强制出来的唯一值。

总结一句：在单相/近单相退化 case 里，g 的最优解集合是一条（甚至一片）“平坦谷”，β 不可辨识；过了 tgt 并不意味着 β 会自动收敛到某个约定值。

如果你希望这类 case 也能“按 β 达标”，就需要额外的规则/正则来打破退化（比如强制小相惩罚、或在评估时把 β 小于某阈值的相视为 inactive 不纳入 e_beta）。你想按哪种方式处理这类退化 case，我可以把 accuracy gate 的规则补进 tp_flash_definitions_and_eval_criteria.tex 并同步实现到 test/tpflash/tp_flash_benchmark.jl。
=#
