"""
TP flash benchmark case definitions.

This file intentionally only defines *cases* (problem instances) and how to compute
their reference solutions (g*, x*, β*). Execution, logging, ranking, and U-score
aggregation live in `test/tpflash/tp_flash_benchmark.jl`.

Cases are derived from `test/test_methods_api_flash.jl` (`@testset "Tp flash algorithms"`),
following `test/tpflash/tp_flash_definitions_and_eval_criteria.tex`.
"""

using Clapeyron

struct TPFlashBenchmarkCase
    id::String
    description::String
    model_builder::Function
    p::Float64
    T::Float64
    feed::Vector{Float64}          # can be normalized z or raw n; benchmark code treats it as "n"
    equilibrium::Symbol            # :auto, :vle, :lle
    numphases::Int                 # optimization assumes a fixed phase count
    logspace::Bool
    g⁺::Union{Float64, Nothing}
    x⁺::Union{Matrix{Float64}, Nothing} # (numphases, numspecies)
    β⁺::Union{Vector{Float64}, Nothing} # length numphases, normalized
end

function tpflash_benchmark_cases()
    cases = TPFlashBenchmarkCase[]

    # --- VLLE (3 phases) ---
    push!(cases, TPFlashBenchmarkCase(
        "vlle_pcsaft_water_cyclohexane_propane_298K_1bar",
        "PCSAFT VLLE (water/cyclohexane/propane), T=298.15K, p=1e5 Pa, z=[0.333,0.333,0.334]",
        () -> PCSAFT(["water", "cyclohexane", "propane"]),
        1e5, 298.15, [0.333, 0.333, 0.334], :auto, 3, false,
        -6.759674475175065,
        [0.9995906526599977 0.00035088514239387127 5.846219760841157e-5; 0.009137394634064861 0.8993697164173008 0.09149288894863436; 0.03333672845431636 0.1224659977187177 0.844197273826966],
        [0.31816821946014534, 0.3210013599067329, 0.36083042063312176],
    ))

    # --- VLE (2 phases) ---
    # NOTE: near single-phase / subcooled degeneracy (from `test/test_methods_api_flash.jl` discussion).
    # Reference has x₁≈x₂≈z and β₂≈0, so g can be near-optimal even when β drifts.
    push!(cases, TPFlashBenchmarkCase(
        "vle_pr_isobutane_nbutane_npenthane_nhexane_282p2K_1bar",
        "PR VLE (IsoButane/n-Butane/n-Pentane/n-Hexane), T=282.2K, p=1e5 Pa, z=[0.25,0.25,0.25,0.25]",
        () -> PR(["IsoButane", "n-Butane", "n-Pentane", "n-Hexane"]),
        1e5, 282.2, [0.25, 0.25, 0.25, 0.25], :vle, 2, false,
        -6.679466011839553,
        [0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25],
        [1.0, 0.0],
    ))

    # NOTE: near bubble/dew boundary: β₂ is small (≈8.9e-4), so the objective can be weakly sensitive to β.
    push!(cases, TPFlashBenchmarkCase(
        "vle_pr_isobutane_nbutane_npenthane_nhexane_282p3K_1bar",
        "PR VLE (IsoButane/n-Butane/n-Pentane/n-Hexane), T=282.3K, p=1e5 Pa, z=[0.25,0.25,0.25,0.25]",
        () -> PR(["IsoButane", "n-Butane", "n-Pentane", "n-Hexane"]),
        1e5, 282.3, [0.25, 0.25, 0.25, 0.25], :vle, 2, false,
        -6.676497175578501,
        [0.24975531522095692 0.2499062625107542 0.25013856689858643 0.2501998553697024; 0.524185521175406 0.3550390728963639 0.09472660107766169 0.02604880485056853],
        [0.9991083897700253, 0.0008916102299747113],
    ))

    # NOTE: near single-phase / saturated-vapor degeneracy (see `test/test_methods_api_flash.jl` logic).
    # Reference has x₁≈x₂≈z and β₁≈0, so g can be near-optimal even when β drifts.
    push!(cases, TPFlashBenchmarkCase(
        "vle_pr_nbutane_npenthane_nhexane_nheptane_450K_1bar",
        "PR VLE (n-butane/n-pentane/n-hexane/n-heptane), T=450K, p=1e5 Pa, z=[0.25,0.25,0.25,0.25]",
        () -> PR(["n-butane", "n-pentane", "n-hexane", "n-heptane"]),
        1e5, 450.0, [0.25, 0.25, 0.25, 0.25], :vle, 2, false,
        -7.279957068200142,
        [0.25 0.25 0.25 0.25; 0.25 0.25 0.25 0.25],
        [0.0, 1.0],
    ))

    push!(cases, TPFlashBenchmarkCase(
        "vle_pr_isobutane_nbutane_npenthane_nhexane_284p4K_1bar",
        "PR VLE (IsoButane/n-Butane/n-Pentane/n-Hexane), T=284.4K, p=1e5 Pa, z=[0.25,0.25,0.25,0.25]",
        () -> PR(["IsoButane", "n-Butane", "n-Pentane", "n-Hexane"]),
        1e5, 284.4, [0.25, 0.25, 0.25, 0.25], :auto, 2, false,
        -6.618441125949686,
        [0.2220766508414057 0.23741130084574857 0.26585578041936375 0.2746562678934819; 0.49736658311503984 0.36152041533284096 0.10953725669926688 0.03157574485285239],
        [0.8985674887273434, 0.10143251127265665],
    ))

    # --- LLE (2 phases) ---
    push!(cases, TPFlashBenchmarkCase(
        "lle_pcsaft_water_cyclohexane_propane_298K_1bar",
        "PCSAFT LLE (water/cyclohexane/propane), T=298.15K, p=1e5 Pa, z=[0.5,0.5,0.0]",
        () -> PCSAFT(["water", "cyclohexane", "propane"]),
        1e5, 298.15, [0.5, 0.5, 0.0], :lle, 2, false,
        -7.577270350886795,
        [0.9996139542909095 0.00038604570909056877 0.0; 0.00970027330568687 0.9902997266943131 0.0],
        [0.49529543445276647, 0.5047045655472335],
    ))

    # --- Activity models (LLE + VLE) ---
    push!(cases, TPFlashBenchmarkCase(
        "lle_unifac_water_hexane_303p15K_1atm",
        "UNIFAC LLE (water/hexane), T=303.15K, p=101325 Pa, z=[0.5,0.5]",
        () -> UNIFAC(["water", "hexane"]),
        101325.0, 303.15, [0.5, 0.5], :lle, 2, false,
        -7.254813021496728,
        [0.999839188463 0.00016081153699992588; 0.009359875304170373 0.9906401246958296],
        [0.49535625648867265, 0.5046437435113273],
    ))

    # NOTE: this case builds a γ-ϕ (GammaPhi) CompositeModel; Clapeyron explicitly
    # disallows VLE evaluation with the `DETPFlash` objective, which is the objective
    # optimized by this benchmark's metaheuristics.
    #=
    push!(cases, TPFlashBenchmarkCase(
        "vle_unifac_octane_heptane_300p15K_2500Pa",
        "UNIFAC VLE (octane/heptane), T=300.15K, p=2500 Pa, z=[0.9,0.1]",
        () -> begin
            if hasfield(UNIFAC, :puremodel)
                UNIFAC(["octane", "heptane"], puremodel = cPR)
            else
                CompositeModel(["octane", "heptane"], liquid = UNIFAC, fluid = cPR)
            end
        end,
        2500.0, 300.15, [0.9, 0.1], :auto, 2, false,
        () -> MichelsenTPFlash(),
    ))
    =#

    # NOTE: correlation-based CompositeModel does not support the DETPFlash objective,
    # so it cannot be benchmarked with the metaheuristic objective used here.
    #=
    push!(cases, TPFlashBenchmarkCase(
        "vle_composite_water_ethanol_358p15K_1atm",
        "CompositeModel VLE (water/ethanol), T=358.15K, p=101325 Pa, z=[0.2,0.8]",
        () -> CompositeModel(["water", "ethanol"], gas = BasicIdeal, liquid = RackettLiquid, saturation = LeeKeslerSat),
        101325.0, 358.15, [0.2, 0.8], :auto, 2, false,
        () -> MichelsenTPFlash(),
        nothing,
        nothing,
        nothing,
    ))
    =#

    # --- Stress / regression cases from #454 (tp_flash2) ---
    function model_454()
        return PR(["n-butane", "n-pentane", "n-hexane", "n-heptane"];
            idealmodel = AlyLeeIdeal,
            userlocations = (;
                Tc = [425.12, 469.7, 507.6, 540.2],
                Pc = [37.96e5, 33.7e5, 30.25e5, 27.4e5],
                Mw = [58.1234, 72.15028, 86.17716, 100.20404],
                acentricfactor = [0.200164, 0.251506, 0.301261, 0.349469],
                k = [
                    0.0        0.0174       -0.0056      0.0033
                    0.0174     0.0          -0.00071726  0.0074
                    -0.0056    -0.00071726   0.0        -0.0078
                    0.0033     0.0074       -0.0078      0.0
                ],
                l = zeros(4, 4),
            ),
        )
    end

    push!(cases, TPFlashBenchmarkCase(
        "vle_454_case1_rr",
        "#454 case 1 (RR): p=153823 Pa, T=321.967K, n=[0.007682,0.9923,1.517e-17,1.918e-31]",
        model_454,
        153_823.0, 321.9670623578307, [0.007682, 0.9923, 1.517e-17, 1.918e-31], :vle, 2, false,
        123.7459904164123,
        [0.00236666244831282 0.9976333375516873 0.0 0.0; 0.007929376297815581 0.9920706237021845 0.0 0.0],
        [0.044445575669629145, 0.9555544243303709],
    ))

    # NOTE: near bubble/dew boundary: β₂ is small (≈1.49e-3), so the objective can be weakly sensitive to β.
    push!(cases, TPFlashBenchmarkCase(
        "vle_454_case2_rr",
        "#454 case 2 (RR): p=701739.83 Pa, T=430.74K, n=[2.984e-14,0.0615,3.48,2.059]",
        model_454,
        701_739.83, 430.74, [2.984e-14, 0.0615, 3.48, 2.059], :vle, 2, false,
        103.70091207385789,
        [5.306960867832837e-15 0.010962897986764673 0.6212133688150436 0.36782373319818634; 1.9503204646204015e-14 0.023231078615878246 0.7284990198696496 0.24826990151445272],
        [0.9985112373218122, 0.0014887626781877376],
    ))

    # NOTE: near single-phase degeneracy (from the #454 block in `test/test_methods_api_flash.jl`):
    # it asserts one phase fraction is exactly zero; under fixed numphases=2 the objective can be flat in β.
    push!(cases, TPFlashBenchmarkCase(
        "vle_454_case3_rr",
        "#454 case 3 (RR): p=1.98555e6 Pa, T=416.663K, n=[55.4614,0.092649,7.265e-9,8.855e-14]",
        model_454,
        1.985550610608908e6, 416.6628781711617, [55.461373286206445, 0.09264900343401582, 7.265116936961075e-9, 8.855321114218425e-14], :vle, 2, false,
        87.73279400794794,
        [0.9983322717803571 0.0016677280888656072 1.3077571409920644e-10 1.5940017928659797e-15; 0.9983322717803571 0.0016677280888656072 1.3077571409920644e-10 1.5940017928659797e-15],
        [0.0, 1.0],
    ))

    # NOTE: near single-phase degeneracy (see `test/test_methods_api_flash.jl` #454 block).
    push!(cases, TPFlashBenchmarkCase(
        "vle_454_case4_rr",
        "#454 case 4 (RR): p=5.35202e5 Pa, T=393.265K, n=[36.4950,0.005799,1.94e-10,2.0e-15]",
        model_454,
        5.35202e5, 393.265, [36.495044786426966, 0.005798955283355085, 1.9416516061189107e-10, 2.0015179988524742e-15], :vle, 2, false,
        92.16253721056026,
        [0.9998411281799807 0.00015887181469993373 5.319470475390174e-12 0.0; 0.9998411281799807 0.00015887181469993373 5.319470475390174e-12 0.0],
        [0.0, 1.0],
    ))

    # NOTE: near single-phase degeneracy (see `test/test_methods_api_flash.jl` #454 block).
    push!(cases, TPFlashBenchmarkCase(
        "vle_454_case5_rr",
        "#454 case 5 (RR): p=442595.319 Pa, T=318.920K, n=[18.6979,9.209e-8,2.317e-22,1.932e-32]",
        model_454,
        442_595.31887270656, 318.91991913774194, [18.697907101753938, 9.208988950434023e-8, 2.317361697667793e-22, 1.9317538045050555e-32], :vle, 2, false,
        114.68417889588606,
        [0.9999999950748557 4.9251442179935655e-9 0.0 0.0; 0.9999999950748557 4.9251442179935655e-9 0.0 0.0],
        [1.0, 0.0],
    ))

    # NOTE: near single-phase degeneracy (see `test/test_methods_api_flash.jl` #454 block).
    push!(cases, TPFlashBenchmarkCase(
        "vle_454_case6_michelsen",
        "#454 case 6 (Michelsen): p=2.20996e6 Pa, T=464.637K, n=[2.756e-6,55.9642,12.8601,1.0820]",
        model_454,
        2.2099578494144413e6, 464.63699168781847, [2.7561794126981888e-6, 55.964211412167195, 12.860133735598001, 1.0819681996211576], :vle, 2, false,
        84.77315106035215,
        [3.9426758071687345e-8 0.8005601572432535 0.18396240071563483 0.015477402614353594; 3.9426758071687345e-8 0.8005601572432535 0.18396240071563483 0.015477402614353594],
        [0.0, 1.0],
    ))

    # NOTE: near single-phase degeneracy (see `test/test_methods_api_flash.jl` #454 block).
    push!(cases, TPFlashBenchmarkCase(
        "vle_454_case7_michelsen",
        "#454 case 7 (Michelsen): p=505777.32 Pa, T=323.960K, n=[75.8306,6.149e-10,0,0]",
        model_454,
        505_777.32016068726, 323.9598978773458, [75.83064821431964, 6.148759359393775e-10, 0.0, 0.0], :vle, 2, false,
        112.95608449327719,
        [0.9999999999918915 8.108541208781074e-12 0.0 0.0; 0.9999999999918915 8.108541208781074e-12 0.0 0.0],
        [1.0, 0.0],
    ))

    return cases
end
