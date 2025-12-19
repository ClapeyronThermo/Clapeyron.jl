## TP flash benchmark (two-stage)

This implements the two-stage workflow described in `test/tpflash/tp_flash_definitions_and_eval_criteria.tex`:

1. **Stage 1** (per algorithm): run one metaheuristic algorithm on all cases and save raw results.
2. **Stage 2** (after all algorithms): load all saved results, compute per-case ranks, U-scores, and aggregate rankings.

### Files

- `test/tpflash/tp_flash_benchmark_cases.jl`: benchmark case definitions (extracted from `test/test_methods_api_flash.jl`).
- In each case, the reference `(g*, x*, β*)` is stored as numeric literals (`g⁺/x⁺/β⁺`) for scoring (no online reference solve).
- `test/tpflash/tp_flash_benchmark_backend.jl`: optimizer backend hook (edit this to swap algorithms by “hacking source”).
- `test/tpflash/tp_flash_benchmark.jl`: evaluation logic (targets, gates, ranking, U-score).
- `test/tpflash/tp_flash_benchmark_stage1.jl`: stage 1 runner.
- `test/tpflash/tp_flash_benchmark_stage2.jl`: stage 2 aggregator.

### Stage 1: run one algorithm and save results

From the repo root:

```bash
julia --project test/tpflash/tp_flash_benchmark_stage1.jl --algo=RDEx --outdir=test/tpflash/results --n=10
```

From the REPL (preferred if you don’t want CLI args):

```julia
include("test/tpflash/tp_flash_benchmark_stage1.jl")
tpflash_benchmark_stage1(; algo="RDEx", outdir="test/tpflash/results", seeds=1:10)
```

Useful options:

- `--algo=<name>`: label used for output filename and for stage 2 comparisons (required).
- `--outdir=<dir>`: where `.ser` result files are written.
- `--n=<int>` / `--seed-start=<int>`: generate seeds `[seed-start, ..., seed-start+n-1]`.
- `--seeds=1,2,3`: explicit seed list.
- `--pop=<int>`: population size for the optimizer backend.
- `--time-limit=<seconds|Inf>`: per-run time limit passed to the optimizer.
- `--stagnation-evals=<int>` and `--stagnation-tol=<float>`: early-stop parameters passed to `RDEx` (0 disables).
- `--filter=<regex>`: only run case ids matching the regex.

This benchmark does not compute reference solutions at runtime; it uses the literals embedded in `tp_flash_benchmark_cases.jl`.

### Swap algorithm (per workflow constraint)

Edit `test/tpflash/tp_flash_benchmark_backend.jl` (or hack the underlying solver implementation) and re-run stage 1 with a new `--algo` label. Repeat for each algorithm variant.

### Stage 2: rank and aggregate after all algorithms finished

```bash
julia --project test/tpflash/tp_flash_benchmark_stage2.jl --indir=test/tpflash/results --outdir=test/tpflash/results
```

From the REPL:

```julia
include("test/tpflash/tp_flash_benchmark_stage2.jl")
tpflash_benchmark_stage2(; indir="test/tpflash/results", outdir="test/tpflash/results")
```

This writes two CSVs into `--outdir`:

- `*_summary.csv`: per-algorithm mean normalized U-score + overall success rate.
- `*_per_case.csv`: per-case U-score and success rate per algorithm.
