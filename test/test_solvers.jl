const SOL = Clapeyron.Solvers

quadratic(x) = x*x - 4
rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
golden_number_fixpoint(x) = one(x) + one(x)/x
@testset "Solvers Module" begin
    @testset "ad_newton" begin
        x0 = 2+rand()
        fdf = SOL.f∂f(quadratic,x0)
        @test fdf[1] ≈ quadratic(x0)
        @test fdf[2] ≈ 2*x0
        @test SOL.ad_newton(quadratic,x0) ≈ 2.0
    end

    @testset "box optimize" begin
        #example of FMinBox in Optim.jl
        lower = [1.25, -2.1]
        upper = [Inf, Inf]
        solution = [1.25, 1.5625]
        initial_x = [2.0, 2.0]
        res = SOL.box_optimize(rosenbrock,initial_x,lower,upper)
        @test all(res.info.minimizer .≈ solution)
        @test res.info.minimum ≈ 0.0625
    end

    @testset "fixpoint" begin
        #example of FMinBox in Optim.jl
        x0 = 2+ rand()
        ϕ = (1+sqrt(5))/2
        @test SOL.fixpoint(golden_number_fixpoint,x0) ≈ ϕ
    end

    function f_diffmcp!(fvec, x)
        fvec[1] = (1-x[1])^2-1.01
    end

    # A difficult MCP.
    #
    # Presented in Miranda and Fackler (2002): "Applied Computational Economics and
    # Finance", p. 51
    #does not converge in NLSolve.jl converges here with NLSolvers.jl
    @testset "nlsolve" begin
        res = SOL.nlsolve(f_diffmcp!,[0.0])
        solution = SOL.x_sol(res)
        v = [1.0]
        f_diffmcp!(v,solution)
        @test v[1] <= sqrt(eps(Float64))
    end
end