@testset "SAFT methods, multi-components" begin
    @printline
    system = PCSAFT(["methanol","cyclohexane"])
    p = 1e5
    T = 313.15
    z = [0.5,0.5]
    p2 = 2e6
    T2 = 443.15
    z2 = [0.27,0.73]

    @testset "Equilibrium properties" begin
        @test Clapeyron.VLLE_pressure(system, T)[1] ≈ 54504.079665621306 rtol = 1E-6
        @test Clapeyron.VLLE_temperature(system, p)[1] ≈ 328.2478837563423 rtol = 1E-6  
    end
    @printline
end
