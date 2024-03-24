using Test
t1 = @elapsed using Clapeyron
using CoolProp #CoolProp ext
using Unitful #Unitful ext
using MultiComponentFlash: MultiComponentFlash

using Clapeyron.LinearAlgebra
using Clapeyron.StaticArrays
using Clapeyron: has_sites,has_groups

#=
This code is copied from ChainRules.jl tests.
I remember that, some time ago, some of the package mantainers said that using JuliaInterpreter in the tests really speed things up.
Clapeyron (compiled) tests take too long, because we are evaluating Each EoS. but in practice, an user will only use a handful of EoS at a time.
=#
@static if Base.VERSION <= v"1.10"
    @eval using JuliaInterpreter #JuliaInterpreter fails on nightly
    union!(JuliaInterpreter.compiled_modules, Any[Base, Base.Broadcast, LinearAlgebra, StaticArrays,Clapeyron.Solvers])
    macro interpret(ex)
        esc(:(JuliaInterpreter.@interpret $ex))
    end
    function include_test(path)
        if isempty(ARGS) || any(occursin(a, path) for a in ARGS)
            println("Testing $path:")  # print so TravisCI doesn't timeout due to no output
            @time Base.include(@__MODULE__(), path) do ex
                Meta.isexpr(ex, :macrocall) && ex.args[1] == Symbol("@testset") || return ex
                return :(@interpret (() -> $ex)())  # interpret testsets using JuliaInterpreter
            end
        else
            # If you provide ARGS like so, then it runs only matching testsets: 
            # Pkg.test("Clapeyron", test_args = ["index", "LinearAlgebra"])
            println("(Not testing $path)")
        end
    end
else
    macro interpret(ex)
        esc(:($ex))
    end
    include_test(path) = include(path)
end



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

@testset "All tests" begin
    include("test_database.jl")
    include("test_solvers.jl")
    include("test_differentials.jl")
    include("test_misc.jl")
    #those two are the main slowdown on the tests.
    include_test("test_models.jl")
    include("test_methods_eos.jl")

    include("test_methods_api.jl")
    include("test_estimation.jl")
    include("test_issues.jl")
end

