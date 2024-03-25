using Test
t1 = @elapsed using Clapeyron
using CoolProp #CoolProp ext
using Unitful #Unitful ext
using MultiComponentFlash: MultiComponentFlash


using Clapeyron.LinearAlgebra
using Clapeyron.StaticArrays
using Clapeyron: has_sites,has_groups

IS_GH = get(ENV,"GITHUB_ACTIONS",false) in ("true", "1", "yes",true,"TRUE")
IS_LOCAL = !IS_GH
IS_STABLE = v"1.6" <= Base.VERSION < v"1.7"
IS_LATEST = v"1.10" <= Base.VERSION < v"1.11"
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
const COVERAGE = !IS_LOCAL && L_NIGHTLY
#ordered by priority
const DISTRIBUTED_WORKER_1 = L_LATEST || W_NIGHTLY
const DISTRIBUTED_WORKER_2 = L_NIGHTLY || A_LATEST
const DISTRIBUTED_WORKER_3 = W_LATEST || A_STABLE
const DISTRIBUTED_WORKER_4 = W_STABLE || A_NIGHTLY
const OTHER_WORKER = !DISTRIBUTED_WORKER_1 && !DISTRIBUTED_WORKER_2 && !DISTRIBUTED_WORKER_3 && !DISTRIBUTED_WORKER_4 && !COVERAGE

DISTRIBUTED_NUMBER = if DISTRIBUTED_WORKER_1
    1
elseif DISTRIBUTED_WORKER_2
    2
elseif DISTRIBUTED_WORKER_3
    3
elseif DISTRIBUTED_WORKER_4
    4
elseif COVERAGE
    -1
else
    0
end
local_str = "Running in " * ifelse(IS_LOCAL,"local","CI") * " mode. " * ifelse(iszero(DISTRIBUTED_NUMBER),"","Distributed worker number: $DISTRIBUTED_NUMBER") * ifelse(DISTRIBUTED_NUMBER == 5," (Coverage)","")

println("""
___________________

Clapeyron.jl tests

$local_str

____________________

""")
#we run coverage in this one:

@info "Loading Clapeyron took $(round(t1,digits = 2)) seconds"
@info "Coolprop: $(Clapeyron.is_coolprop_loaded())"
#Disable showing citations
ENV["CLAPEYRON_SHOW_REFERENCES"] = "FALSE"

macro printline()  # useful in hunting for where tests get stuck
    file = split(string(__source__.file), "/")[end]
    printstyled(">>", file, ":", __source__.line, "\n", color=:light_black)
end

#fix to current tests
function GERG2008(components;verbose = false,reference_state = nothing)
    return MultiFluid(components;
    mixing = AsymmetricMixing,
    departure = EmpiricDeparture,
    pure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/pures"],
    mixing_userlocations  = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/mixing/GERG2008_mixing_unlike.csv"],
    departure_userlocations = String["@REMOVEDEFAULTS","@DB/Empiric/GERG2008/departure/GERG2008_departure_unlike.csv"],
    reference_state = reference_state,
    coolprop_userlocations = false,
    verbose = verbose,
    Rgas = Clapeyron.R̄)
end

function test_gibbs_duhem(model,V,T,z;rtol = 1e-14)
    _,G,∑μᵢzᵢ = Clapeyron.gibbs_duhem(model,V,T,z)
    @test G ≈ ∑μᵢzᵢ rtol = rtol
end

if DISTRIBUTED_WORKER_4 || IS_LOCAL || COVERAGE || OTHER_WORKER
    include("test_database.jl")
    include("test_solvers.jl")
    include("test_differentials.jl")
    include("test_misc.jl")
    include("test_models_saft_pc.jl")
end

if DISTRIBUTED_WORKER_3 || IS_LOCAL || COVERAGE || OTHER_WORKER
    include("test_models_cubic.jl")
    include("test_models_saft_others.jl")
end

if DISTRIBUTED_WORKER_2 || IS_LOCAL || COVERAGE || OTHER_WORKER
    include("test_models_others.jl")
end
if DISTRIBUTED_WORKER_1 || IS_LOCAL || COVERAGE || OTHER_WORKER
    include("test_models_saft_vr.jl")
end
if DISTRIBUTED_WORKER_4 || IS_LOCAL || COVERAGE || OTHER_WORKER
    include("test_methods_eos.jl")
end
if DISTRIBUTED_WORKER_3 || IS_LOCAL || COVERAGE || OTHER_WORKER
    include("test_methods_api.jl")
end
if DISTRIBUTED_WORKER_2 || IS_LOCAL || COVERAGE || OTHER_WORKER
    include("test_estimation.jl")
end

if DISTRIBUTED_WORKER_1 || IS_LOCAL || COVERAGE || OTHER_WORKER
    include("test_issues.jl")
end
