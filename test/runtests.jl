using Test, Unitful

t1 = @elapsed using Clapeyron
@info "Loading Clapeyron took $(round(t1,digits = 2)) seconds"

#Disable showing citations
ENV["CLAPEYRON_SHOW_REFERENCES"] = "FALSE"

macro printline()  # useful in hunting for where tests get stuck
    file = split(string(__source__.file), "/")[end]
    printstyled("  ", file, ":", __source__.line, "\n", color=:light_black)
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
end
