using OpenSAFT, Test

@testset "All tests" begin
    include("test_database.jl")
    include("test_models.jl")
    include("test_methods.jl")
end
