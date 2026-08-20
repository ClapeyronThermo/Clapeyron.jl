const SOL = Clapeyron.Solvers
@printline
quadratic(x) = x*x - 4
rosenbrock(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
golden_number_fixpoint(x) = one(x) + one(x)/x
fgh_lgmx(x) = (log(x)+x,1/x +1,-1/(x^2))
fg_2(x,y) = (x+2*y)^2
function quadratic_fixpoint(y,x)
    y[1] = sqrt(-x[2]^2 + 5)
    y[2] = 2*sqrt(-x[1]^2 + 9)/3
    y
end

minlog(x) = ((x - 2)^2 - log(x))

#=
2*(x - 2) - 1/x = 0
2x - 4 - 1/x = 0
2x^2 - 4x - 1 = 0
x^2 - 2x - 0.5 = 0
(2 + sqrt(4 + 2))/2
(2 + sqrt(6))/2
=#

function opt_test(k,method)
    my_rosenbrock(x) = (1.0 - k*x[1])^2 + 100 * (x[2] - x[1]^2)^2
    x0 = Clapeyron.SA[0.0,0.0]
    return SOL.solution(SOL.optimize(my_rosenbrock,x0,method)) #uses newton as default.
    #return solution(optimize(rosenbrock,x0,LineSearch(Newton()),OptimizationOptions())
end


function opt_test(k)
    my_rosenbrock(x) = (1.0 - k*x[1])^2 + 100 * (x[2] - x[1]^2)^2
    x0 = Clapeyron.SA[0.0,0.0]
    return SOL.solution(SOL.optimize(my_rosenbrock,x0,SOL.LineSearch(SOL.Newton()))) #uses newton as default.
    #return solution(optimize(rosenbrock,x0,LineSearch(Newton()),OptimizationOptions())
end
@testset verbose = true "Solvers Module" begin
    @testset "Newton, Halley" begin
        x0 = 2+rand()
        fdf = SOL.f∂f(quadratic,x0)
        @test fdf[1] ≈ quadratic(x0)
        @test fdf[2] ≈ 2*x0
        @test SOL.ad_newton(quadratic,x0) ≈ 2.0
        @test SOL.newton(x->(quadratic(x),2*x),x0) ≈ 2.0
        @test SOL.halley(fgh_lgmx,0.5)≈ 0.567143290409784
    end

    @testset "fixpoint" begin
        #example of FMinBox in Optim.jl
        x0 = 2+ rand()
        ϕ = (1+sqrt(5))/2
        x1 = [0.1,0.1]
        @test @inferred(SOL.fixpoint(golden_number_fixpoint,x0)) ≈ ϕ
        @test @inferred(SOL.fixpoint(golden_number_fixpoint,x0,SOL.AitkenFixPoint())) ≈ ϕ
        @test @inferred(SOL.fixpoint(quadratic_fixpoint,x1)) ≈ [0.6*sqrt(5),0.8*sqrt(5)]
    end

    function f_diffmcp!(fvec, x)
        fvec[1] = (1-x[1])^2-1.01
    end

    @testset "roots3" begin
        #@test [@inferred(SOL.roots3(SA[-6im, -(3 + 4im), 2im-2, 1.0]))...] ≈ [3, -2im, -1]
        @test @inferred(SOL.roots3(1.0, -3.0, 3.0, -1.0)) ≈ [1, 1, 1]
        @test @inferred(SOL.roots3([1.0, -3.0, 3.0, -1.0])) ≈ [1, 1, 1]
        
        p = NTuple{4,Float64}[]
        s = Tuple{Int64,Float64,Float64,Float64}[]
        push!(p, (0.4179831841886642, -16.86256511617483, 1.6508521434627874, 1.0))
        push!(p, (0.7336510424849684, -17.464843881306653, 1.6925644348171853, 1.0))
        push!(p, (-1, 3, -3, 1))
        push!(p, (1, -1, -1, 1))
        push!(p, (5.025484622448689e22, -1.8479625969997355e28, 1.6619181073681473e33, -4.259242533341054e37))
        
        push!(s, (3, -5.0238891763914015, 0.024849000208331546, 3.3481880327202824))
        push!(s, (3, -5.126952070514365, 0.04218405979418538, 3.392203575902995))
        push!(s, (1, 1.0, 1.0, 1.0))
        push!(s, (2, -1.0, 1.0, 1.0))
        push!(s, (3, 4.02794612966951e-6, 1.3867007835810342e-5, 2.1124146129402466e-5))

        for i in 1:length(p)
            sol = @inferred(SOL.real_roots3(p[i]))
            soli = s[i]
            @test sol[1] == soli[1]
            v = [sol[2],sol[3],sol[4]]
            vi = [soli[2],soli[3],soli[4]]
            @test v ≈ vi
        end
        @test collect(SOL.real_roots3(1, -1, -1, 1)) ≈ [2, -1.0, 1.0, 1.0] # one doubly degenerate
        @test SOL.real_roots3(1, 0, 1, 2) == (1, -1.0, -1.0, -1.0) # only one real root
    end

    # A difficult MCP.
    #
    # Presented in Miranda and Fackler (2002): "Applied Computational Economics and
    # Finance", p. 51
    #does not converge in NLSolve.jl converges here with NLSolvers.jl
    @testset "nlsolve" begin
        res = SOL.nlsolve(f_diffmcp!,[0.0],SOL.NLSolvers.TrustRegion(SOL.NLSolvers.Newton()),SOL.NLSolvers.NEqOptions(),ForwardDiff.Chunk{1}())
        res2 = SOL.nlsolve(f_diffmcp!,[0.0],SOL.NLSolvers.LineSearch(SOL.NLSolvers.Newton()),SOL.NLSolvers.NEqOptions(),ForwardDiff.Chunk{1}())
        @test SOL.x_sol(res) isa Vector
        @test SOL.x_sol(res2) isa Vector
        solution = SOL.x_sol(res)
        v = [1.0]
        f_diffmcp!(v,solution)
        @test v[1] <= sqrt(eps(Float64))
    end

    @testset "det_22" begin
        a1 = 1
        a2 = 2
        a3 = 3
        a4 = 4
        @test SOL.det_22(1,2,3,4) == a1*a2 - a3*a4
    end

    @testset "AD - misc" begin
        @test SOL.gradient2(fg_2,1,1) == [6,12]
        @test eltype(SOL.gradient2(fg_2,1.0,1)) == Float64
        @test eltype.(SOL.∂2(fg_2,1.0,1)) == (Float64, Float64, Float64)
        ad = SOL.ADScalarObjective(rosenbrock,zeros(2))
        @test ad.f(zeros(2)) == 1.0
        @test ad.fg(ones(2),zeros(2))[2] == [-2.0,0.0]
        @test ad.fgh(ones(2),ones(2,2),zeros(2))[3] == [2.0 0.0; 0.0 200.0]

        d1 = SOL.grad_at_i(rosenbrock,[2.0,2.0],1)
        d2 = SOL.grad_at_i(rosenbrock,[2.0,2.0],2)
        d = SOL.gradient(rosenbrock,[2.0,2.0])
        @test [d1,d2] == d
        
        fx = SOL.strong_zero(Clapeyron.ForwardDiff.Dual(1.2,0.5)) do x
            sqrt(x)
        end
        @test Clapeyron.primalval(fx) ≈ sqrt(1.2)
    end

    @testset "evalexppoly" begin
        @test Clapeyron.evalexppoly(2,(1,2,3),(3,2,1)) == 22
    end

    @testset "optimize" begin
        xs = [1/3,1/9] #solution
        x1 = opt_test(3.0)
        @test x1 ≈ xs

        x2 = opt_test(3.0,SOL.NelderMead())
        @test x2 ≈ xs rtol = 1e-3

        xs_1v = (2 + sqrt(6))/2
        x3 = SOL.solution(SOL.optimize(minlog,(1.5,2.5),SOL.BrentMin()))
        @test x3 ≈ xs_1v

        x4 = SOL.solution(SOL.optimize(minlog,(1.5,2.5),SOL.BoundOptim1Var()))
        @test x4 ≈ xs_1v
    end

end
@printline
