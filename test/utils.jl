IS_GH = get(ENV,"GITHUB_ACTIONS",false) in ("true", "1", "yes",true,"TRUE")
IS_LOCAL = !IS_GH
IS_STABLE = v"1.10" <= Base.VERSION < v"1.11"
IS_LATEST = v"1.12" <= Base.VERSION < v"1.12.9999"
IS_OTHER = !IS_STABLE && !IS_LATEST
const W_STABLE = IS_STABLE && (Base.Sys.iswindows())
const W_LATEST = IS_LATEST && (Base.Sys.iswindows())
const W_NIGHTLY = IS_OTHER && (Base.Sys.iswindows())
const L_STABLE = IS_STABLE && (Base.Sys.islinux())
const L_LATEST = IS_LATEST && (Base.Sys.islinux())
const L_NIGHTLY = IS_OTHER && (Base.Sys.islinux())
const A_STABLE = IS_STABLE && (Base.Sys.isapple())
const A_LATEST = IS_LATEST && (Base.Sys.isapple())
const A_NIGHTLY = IS_OTHER && (Base.Sys.isapple())
#ordered by priority
const DISTRIBUTED_WORKER_1 = L_NIGHTLY || W_LATEST
const DISTRIBUTED_WORKER_2 = L_STABLE || A_LATEST
const DISTRIBUTED_WORKER_3 = W_NIGHTLY || A_STABLE
const DISTRIBUTED_WORKER_4 = W_STABLE || A_NIGHTLY
const DISTRIBUTED_WORKER_5 = L_LATEST

const OTHER_WORKER = !DISTRIBUTED_WORKER_1 && !DISTRIBUTED_WORKER_2 && !DISTRIBUTED_WORKER_3 && !DISTRIBUTED_WORKER_4 && !DISTRIBUTED_WORKER_5

DISTRIBUTED_NUMBER = if DISTRIBUTED_WORKER_1
    1
elseif DISTRIBUTED_WORKER_2
    2
elseif DISTRIBUTED_WORKER_3
    3
elseif DISTRIBUTED_WORKER_4
    4
elseif DISTRIBUTED_WORKER_5
    5
else
    0
end

local_str = "Running in " * ifelse(IS_LOCAL,"local","CI") * " mode. " * ifelse(iszero(DISTRIBUTED_NUMBER),"","Distributed worker number: $DISTRIBUTED_NUMBER")

println("""
___________________

Clapeyron.jl tests

$local_str
Running all tests in all workers = $ALL_TESTS
____________________

""")

function include_distributed(path,value::Int)
    workers = (DISTRIBUTED_WORKER_1,DISTRIBUTED_WORKER_2,DISTRIBUTED_WORKER_3,DISTRIBUTED_WORKER_4,DISTRIBUTED_WORKER_5)
    if value in 1:5
        worker = workers[value]
    else
        worker = false
    end
    if worker || IS_LOCAL || OTHER_WORKER || ALL_TESTS
        Base.include(@__MODULE__(), path)
    end
end

macro printline()  # useful in hunting for where tests get stuck
    file = split(string(__source__.file), Base.Filesystem.path_separator)[end]
    printstyled(">>", file, ":", __source__.line, "\n", color=:light_black)
end

function model500()

#created via eos_repr
return SAFTgammaMie{BasicIdeal, Float64, Int64}(
        ["DIMETHYLAMINE"], 
        GroupParam(["DIMETHYLAMINE"], [["NH", "CH3"]], :unknown, [[1, 2]], [Matrix{Int64}(undef, 0, 0)], [[2, 1]], ["CH3", "NH"], [[2, 1]], String[]), 
        SiteParam(["CH3", "NH"], [String[], ["H", "e1"]], Clapeyron.PackedVofV([1, 1, 3], [1, 1]), [Int64[], [1, 2]], ["H", "e1", "e2"], [[0, 0, 0], [1, 1, 0]], [[0, 0, 0], [1, 2, 0]], String[], nothing), 
        Clapeyron.SAFTgammaMieParam{Float64}(
        SingleParam{Int64}("vst", ["CH3", "NH"], [1, 1], Bool[0, 0], ["/home/22796002/PhD/Smarts/SAFTgammaMie_like_Clapeyron.csv"], String[]), 
        SingleParam{Float64}("S", ["CH3", "NH"], [0.31696615176726495, 0.15785492572578758], Bool[0, 0], nothing, nothing), 
        PairParam{Float64}("lambda_a", ["CH3", "NH"], [6.0 6.0; 6.0 6.0], Bool[0 1; 1 0], nothing, nothing), 
        PairParam{Float64}("lambda_r", ["CH3", "NH"], [31.239630776055115 8.307492753603125; 8.307492753603125 18.449975142848835], Bool[0 0; 0 0], nothing, nothing), 
        PairParam{Float64}("sigma", ["CH3", "NH"], [5.245968421264455e-10 4.0342843067136425e-10; 4.0342843067136425e-10 2.822600192162831e-10], Bool[0 1; 1 0], nothing, nothing), 
        PairParam{Float64}("epsilon", ["CH3", "NH"], [467.0612223874854 268.98688435945206; 268.98688435945206 374.46738001552916], Bool[0 0; 0 0], nothing, nothing), 
        AssocParam{Float64}("epsilon_assoc", ["CH3", "NH"], 
        Clapeyron.Compressed4DMatrix{Float64, Vector{Float64}}([1693.6574275348516], [(2, 2)], [(1, 2)], (2, 2), (3, 3)), [String[], ["H", "e1"]], nothing, nothing), 
        AssocParam{Float64}("bondvol", ["CH3", "NH"], 
        Clapeyron.Compressed4DMatrix{Float64, Vector{Float64}}([4.532862604822759e-29], [(2, 2)], [(1, 2)], (2, 2), (3, 3)), [String[], ["H", "e1"]], nothing, nothing), 
        MixedGCSegmentParam{Float64}("mixed segment", ["DIMETHYLAMINE"], Clapeyron.PackedVofV([1, 3], [0.6339323035345299, 0.15785492572578758]))), 
        BasicIdeal(), 
        SAFTVRMie{BasicIdeal, Float64}(["DIMETHYLAMINE"], 
        SiteParam(["DIMETHYLAMINE"], [["NH/H", "NH/e1"]], Clapeyron.PackedVofV([1, 3], [1, 1]), [[21, 22]], ["CH3OH/H", "CH3OH/e1", "CO2/e1", "CO2/e2", "COOH/H", "COOH/e1", "COOH/e2", "CH2OH/H", "CH2OH/e1", "COO/e1", "OH/H", "OH/e1", "cNH/H", "cNH/e1", "NH2/H", "NH2/e1", "CH3CO/e1", "CHOH/H", "CHOH/e1", "N/e1", "NH/H", "NH/e1", "cN/e1", "aCCOOH/H", "aCCOOH/e1", "aCCOOH/e2", "aCH/e1", "aCOH/H", "aCOH/e1", "aCOH/e2", "aCCH3/e1", "aCCH2/e1", "H2O/H", "H2O/e1", "CH3COCH3/H", "CH3COCH3/e1", "CH3COCH3/e2"], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], String[], [[(23, 1), (23, 2)]]), 
        Clapeyron.SAFTVRMieParam{Float64}(
        SingleParam{Float64}("Mw", ["DIMETHYLAMINE"], [45.085], Bool[0], nothing, nothing), 
        SingleParam{Float64}("segment", ["DIMETHYLAMINE"], [0.7917872292603174], Bool[0], nothing, nothing), 
        PairParam{Float64}("sigma", ["DIMETHYLAMINE"], [4.854448918280263e-10;;], Bool[0;;], nothing, nothing), 
        PairParam{Float64}("lambda_a", ["DIMETHYLAMINE"], [6.000000000000001;;], Bool[0;;], nothing, nothing), 
        PairParam{Float64}("lambda_r", ["DIMETHYLAMINE"], [23.41048565515027;;], Bool[0;;], nothing, nothing), 
        PairParam{Float64}("epsilon", ["DIMETHYLAMINE"], [400.14816241181626;;], Bool[0;;], nothing, nothing), 
        AssocParam{Float64}("epsilon_assoc", ["DIMETHYLAMINE"], 
        Clapeyron.Compressed4DMatrix{Float64, Vector{Float64}}([1693.6574275348516], [(1, 1)], [(1, 2)], (1, 1), (7, 7)), [["NH/H", "NH/e1"]], nothing, nothing), 
        AssocParam{Float64}("bondvol", ["DIMETHYLAMINE"], 
        Clapeyron.Compressed4DMatrix{Float64, Vector{Float64}}([4.532862604822759e-29], [(1, 1)], [(1, 2)], (1, 1), (7, 7)), [["NH/H", "NH/e1"]], nothing, nothing)), 
        BasicIdeal(), 
        AssocOptions(1.0e-12, 1.0e-12, 1000, 0.5, :nocombining), ["10.1063/1.4819786", "10.1080/00268976.2015.1029027"]), :default, 
        AssocOptions(1.0e-12, 1.0e-12, 1000, 0.5, :nocombining), ["10.1063/1.4851455", "10.1021/je500248h"])
end
