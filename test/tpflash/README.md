## TP flash benchmark (two-stage)

This folder implements the two-stage workflow described in the benchmark documentation attached to [PR #516](https://github.com/ClapeyronThermo/Clapeyron.jl/pull/516).

1. **Stage 1** (per algorithm): run one optimizer on all cases and save raw run logs (`.ser`).
2. **Stage 2** (after all algorithms): load all logs, compute per-case ranks and U-scores, and write CSV summaries.

### Entry point (recommended)

Use `test/tpflash/run_benchmark_suites.jl` as the main entry point. It runs a tuning suite (population-size sweep) and a validation suite, and then triggers stage 2 aggregation.

From the repo root:

```bash
julia --project test/tpflash/run_benchmark_suites.jl
```

From the REPL:

```julia
include("test/tpflash/run_benchmark_suites.jl")
```

Edit `run_benchmark_suites.jl` to adjust:
- seed ranges (tuning vs validation)
- population sizes and `algo` labels (e.g. `sass15`)
- FE budget (`fe_factor`) and target tolerances (`atol_tgt`, `rtol_tgt`)

Note: `backend=:bbo` requires `BlackBoxOptim.jl` available in your active environment. If you do not have it, comment out the `:bbo` suite in `run_benchmark_suites.jl`.

### Files

- `test/tpflash/tp_flash_benchmark_cases.jl`: case definitions extracted from `test/test_methods_api_flash.jl`. Each case stores reference literals `(g*, x*, beta*)` for scoring (no online reference solve).
- `test/tpflash/tp_flash_benchmark_backend.jl`: optimizer backend hook (edit this to map `backend` to a concrete optimizer).
- `test/tpflash/tp_flash_benchmark.jl`: scoring logic (targets, gates, multi-tag failures, ranking, U-score aggregation).
- `test/tpflash/tp_flash_benchmark_stage1.jl`: stage 1 runner (writes one `.ser` per algorithm label).
- `test/tpflash/tp_flash_benchmark_stage2.jl`: stage 2 aggregator (reads `.ser`, writes CSV summaries).
- `test/tpflash/run_benchmark_suites.jl`: convenience driver that runs stage 1 suites + stage 2.

### Stage 1 (manual)

From the repo root:

```bash
julia --project test/tpflash/tp_flash_benchmark_stage1.jl --algo=sass15 --backend=sass --outdir=test/tpflash/results --n=10
```

From the REPL:

```julia
include("test/tpflash/tp_flash_benchmark_stage1.jl")
tpflash_benchmark_stage1(; algo="sass15", backend=:sass, outdir="test/tpflash/results", seeds=1:10)
```

Common options:
- `--algo=<label>`: label used for output filename and stage 2 comparisons (required).
- `--backend=<sass|bbo>`: optimizer backend selector.
- `--outdir=<dir>`: where `.ser` logs are written.
- `--n=<int>` / `--seed-start=<int>`: generate seeds `[seed-start, ..., seed-start+n-1]`.
- `--seeds=1,2,3`: explicit seed list.
- `--pop=<int>`: population size.
- `--fe-factor=<int>`: FE budget scaling factor (`FE_max = fe_factor * D`).
- `--time-limit=<seconds|Inf>`: per-run time limit.
- `--stagnation-evals=<int>` and `--stagnation-tol=<float>`: backend stagnation stop controls (0 disables).
- `--filter=<regex>`: only run case ids matching the regex.

### Stage 2

```bash
julia --project test/tpflash/tp_flash_benchmark_stage2.jl --indir=test/tpflash/results --outdir=test/tpflash/results
```

This writes two CSVs into `--outdir`:
- `*_summary.csv`: per-algorithm mean normalized U-score and overall success rate.
- `*_per_case.csv`: per-case U-score and per-case success rate per algorithm.
