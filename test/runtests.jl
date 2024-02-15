using Test
t1 = @elapsed using Clapeyron
using CoolProp #CoolProp ext
using Unitful #Unitful ext
using MultiComponentFlash: MultiComponentFlash

@info "Loading Clapeyron took $(round(t1,digits = 2)) seconds"
@info "Coolprop: $(Clapeyron.is_coolprop_loaded())"
#Disable showing citations
ENV["CLAPEYRON_SHOW_REFERENCES"] = "FALSE"

macro printline()  # useful in hunting for where tests get stuck
    file = split(string(__source__.file), "/")[end]
    printstyled("  ", file, ":", __source__.line, "\n", color=:light_black)
end

#fix to current tests
function GERG2008(components::Vector{String};verbose = false,reference_state = nothing)
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
    include("test_models.jl")
    include("test_methods_eos.jl")
    include("test_methods_api.jl")
    include("test_estimation.jl")
    include("test_issues.jl")
end
