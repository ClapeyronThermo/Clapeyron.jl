include("tp_flash_benchmark_stage1.jl");
include("tp_flash_benchmark_stage2.jl");
seeds = -20:-1;
fe_factor = 3000;
time_limit = Inf;
atol_tgt, rtol_tgt = 1e-12, 1e-12;

## Tuning
for pop in [10, 15, 20, 25, 30, 40, 50, 60, 70, 80]
    tpflash_benchmark_stage1(; algo="bbo$pop", backend=:bbo, seeds, population_size=pop, fe_factor, time_limit, stagnation_evals=0, stagnation_tol=0.0, filter=nothing, atol_tgt, rtol_tgt)
end
for pop in [10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 80]
    tpflash_benchmark_stage1(; algo="sass$pop", backend=:sass, seeds, population_size=pop, fe_factor, time_limit, stagnation_evals=0, stagnation_tol=0.0, filter=nothing, atol_tgt, rtol_tgt)
end

## Validation
seeds = 1:20;
pop = 20;
tpflash_benchmark_stage1(; algo="bbo$pop", backend=:bbo, seeds, population_size=pop, fe_factor, time_limit, stagnation_evals=0, stagnation_tol=0.0, filter=nothing, atol_tgt, rtol_tgt)
pop = 15;
tpflash_benchmark_stage1(; algo="sass$pop", backend=:sass, seeds, population_size=pop, fe_factor, time_limit, stagnation_evals=0, stagnation_tol=0.0, filter=nothing, atol_tgt, rtol_tgt)

##----
tpflash_benchmark_stage2(indir="test/tpflash/results/tuning/sass")
tpflash_benchmark_stage2(indir="test/tpflash/results/validation")
