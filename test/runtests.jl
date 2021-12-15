using Clapeyron, Test, Unitful

@testset "All tests" begin
    include("test_database.jl")
    include("test_solvers.jl")
    include("test_differentials.jl")
    include("test_misc.jl")
    include("test_models.jl")
    include("test_methods.jl")
end
