include("tp_flash_benchmark_stage1.jl");
include("tp_flash_benchmark_stage2.jl");
include("plot_tpflash_summary.jl");
seeds = -20:-1;
fe_factor = 3000;
time_limit = Inf;
atol_tgt, rtol_tgt = 1e-12, 1e-12;

##----
clear()
for pop in [20]
    tpflash_benchmark_stage1(; algo="bbo$pop", backend=:bbo, seeds, population_size=pop, fe_factor, time_limit, stagnation_evals=0, stagnation_tol=0.0, filter=nothing, atol_tgt, rtol_tgt)
end
for pop = [40]
    tpflash_benchmark_stage1(; algo="rdex$pop", backend=:rdex, seeds, population_size=pop, fe_factor, time_limit, stagnation_evals=0, stagnation_tol=0.0, filter=nothing, atol_tgt, rtol_tgt)
end
for pop = [10]
    tpflash_benchmark_stage1(; algo="sass$pop", backend=:sass, seeds, population_size=pop, fe_factor, time_limit, stagnation_evals=0, stagnation_tol=0.0, filter=nothing, atol_tgt, rtol_tgt)
end

##----
tpflash_benchmark_stage2(indir="test/tpflash/results/tuning/sass")
tpflash_benchmark_stage2(indir="test/tpflash/results/validation")

##----
path = "D:\\Julia-packages\\Clapeyron.jl\\test\\tpflash\\results\\tuning\\rdex\\"
plot_tpflash_summary_csv(path * "tpflash_benchmark_2025-12-19_140721_summary.csv",
    outpath=path * "plot.pdf", title="RDEx turning")
