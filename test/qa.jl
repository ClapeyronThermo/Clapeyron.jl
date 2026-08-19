#File to add any quality assurance tests.

@printline

@testset verbose = true "QA" begin
    @testset "ambiguities" begin
        ambiguities = Test.detect_ambiguities(Clapeyron)
        if !isempty(ambiguities)
            foreach(display, ambiguities)
        end
        @test length(ambiguities) == 0
    end
end
