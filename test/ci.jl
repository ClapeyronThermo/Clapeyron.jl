IS_GH = get(ENV, "GITHUB_ACTIONS", false) in ("true", "1", "yes", true, "TRUE")
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

const DISTRIBUTED_NUMBER = if DISTRIBUTED_WORKER_1
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

"""
    CI_ALL_TESTS

if set to `true`, run all tests on all workers in CI.
"""
CI_ALL_TESTS = false

function __run_test(value::Int)
    return (DISTRIBUTED_NUMBER == value) || IS_LOCAL || OTHER_WORKER || CI_ALL_TESTS
end

local_str = "Running in " * ifelse(IS_LOCAL, "local", "CI") * " mode. " * ifelse(iszero(DISTRIBUTED_NUMBER), "", "Distributed worker number: $DISTRIBUTED_NUMBER")

println("""
___________________

Clapeyron.jl tests

$local_str
Running all tests in all workers = $CI_ALL_TESTS
____________________

""")
