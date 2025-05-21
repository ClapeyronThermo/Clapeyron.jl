using Clapeyron, Test, NLSolvers, ForwardDiff

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


@testset "Solvers Module" begin
    @testset "newton,halley" begin
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
        @test @inferred(SOL.real_roots3((0.4179831841886642, -16.86256511617483, 1.6508521434627874, 1.0))) ==(3, -5.0238891763914015, 0.02484900020833184, 3.3481880327202824)
        @test @inferred(SOL.real_roots3(0.7336510424849684, -17.464843881306653, 1.6925644348171853, 1.0)) == (3, -5.126952070514365, 0.04218405979418567, 3.392203575902995)
        @test @inferred(SOL.real_roots3(-1, 3, -3, 1)) == (2, 1.0, 1.0, 1.0) # triply degenerate root
        @test collect(SOL.real_roots3(1, -1, -1, 1)) ≈ [2, -1.0, 1.0, 1.0] # one doubly degenerate
        @test SOL.real_roots3(1, 0, 1, 2) == (1, -1.0, -1.0, -1.0) # only one real root
    end

    # A difficult MCP.
    #
    # Presented in Miranda and Fackler (2002): "Applied Computational Economics and
    # Finance", p. 51
    #does not converge in NLSolve.jl converges here with NLSolvers.jl
    @testset "nlsolve" begin
        res = SOL.nlsolve(f_diffmcp!,[0.0],TrustRegion(Newton()),SOL.NLSolvers.NEqOptions(),ForwardDiff.Chunk{1}())
        res2 = SOL.nlsolve(f_diffmcp!,[0.0],LineSearch(Newton()),SOL.NLSolvers.NEqOptions(),ForwardDiff.Chunk{1}())
        @test SOL.x_sol(res) isa Vector
        @test SOL.x_sol(res2) isa Vector
        solution = SOL.x_sol(res)
        v = [1.0]
        f_diffmcp!(v,solution)
        @test v[1] <= sqrt(eps(Float64))
    end


    @testset "Rachford Rice" begin
    #error, all Ks > 0, from CalebBell/Chemicals
    zs = [0.0, 0.0885053990596404, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.721469918219507, 0.01742948033685831,
            0.1725952023839942]
    Ks = [192.3625321718047, 105.20070698573475, 76.30532397111489, 37.090890262982, 21.38862102676539, 18.093547012968767,
            10.319129068837443, 9.001962137200403, 4.565198737490148, 340.69153314749224, 269.09234343328467,
            21.858385052861507]
    @test Clapeyron.rr_vle_vapor_fraction(Ks,zs) == Inf

    #here are some tests, from the paper.
    #the paper has errors, z2 has 10 components, 9 and 10 repeated.
    #and after that, they dont include the bmax,bmin in their selected pair.
    z1 = [0.30,0.15, 0.05, 0.02, 0.01, 0.02, 0.02, 0.03, 0.07, 0.33]
    k1 = [3.0E+00, 2.0E+00, 1.1E+00, 8.0E-01, 4.0E-01, 1.0E-01, 5.0E-02, 2.0E-02, 1.0E-02, 1.0E-04]
    idxs1,ks1,zs1 = Clapeyron.rr_find_strongest(k1,z1)
    @test all(in.(idxs1,Ref((1,2,9,10))))
    @test Clapeyron.rr_vle_vapor_fraction(k1,z1) ≈ 0.1769 atol = 1e-4
    z2 = [0.00034825,0.01376300,0.13084000,0.10925000,0.00001000,0.51009000,0.23564000,0.00006000]
    k2 = [5.2663E+02, 5.0400E+01, 1.6463E+00, 8.7450E-01, 1.5589E-01, 3.6588E-02, 2.6625E-02, 4.8918E-06]
    idxs2,ks2,zs2 = Clapeyron.rr_find_strongest(k2,z2)
    #the paper is wrong here, it stablishes (1,2,6,7) ,does not include bmax
    #we still converge to the correct root
    @test all(in.(idxs2,Ref((1,2,6,8))))
    @test Clapeyron.rr_vle_vapor_fraction(k2,z2) ≈ 0.00327 atol = 1e-4

    z3 = [0.0187002, 0.0243002, 0.5419054, 0.0999010, 0.0969010, 0.0400004, 0.0212002, 0.0148001, 0.0741507, 0.0350404, 0.0173602, 0.0157402]
    k3 = [1.32420, 1.12778, 1.22222, 1.11760, 9.88047E-01, 8.94344E-01, 7.87440E-01, 7.43687E-01, 8.11797E-01, 6.93279E-01, 5.09443E-01, 2.28721E-01]
    idxs3,ks3,zs3 = Clapeyron.rr_find_strongest(k3,z3)
    @test all(in.(idxs3,Ref((1,3,9,12))))
    @test Clapeyron.rr_vle_vapor_fraction(k3,z3) ≈ 0.98725 atol = 1e-4

    #testing exact solver,3 comps:
    zs_3 = [0.5, 0.3, 0.2]
    Ks_3 = [1.685, 0.742, 0.532]
    xs_3 = [0.33940869696634357, 0.3650560590371706, 0.2955352439964858]
    ys_3 = [0.5719036543882889, 0.27087159580558057, 0.15722474980613044]
    β_3 = 0.6907302627738544
    β_3_sol = Clapeyron.rr_vle_vapor_fraction(Ks_3,zs_3)
    @test  β_3_sol ≈ β_3
    @test Clapeyron.rr_flash_eval(Ks_3,zs_3,β_3_sol) <= 4*eps(β_3)
    @test Clapeyron.rr_flash_vapor(Ks_3,zs_3,β_3_sol) ≈ ys_3
    @test Clapeyron.rr_flash_liquid(Ks_3,zs_3,β_3_sol) ≈ xs_3

    #Extreme Kvalues, it should give β = 0.99999
    Ks_extreme = [15.464909530837806, 7006.64008090944, 1.8085837711444488, 0.007750676421035811, 30.98450366497431]
    zs_extreme = [0.26562380186293233, 0.04910234829974003, 0.284394553603828, 0.006300023876072552, 0.3945792723574272]
    @test Clapeyron.rr_vle_vapor_fraction(Ks_extreme,zs_extreme) ≈ 0.999999999

    # Case where the evaluated point is right on the boundary
    Ks_eps_0 = [1.2566703532018493e-21, 3.3506275205339295, 1.0300675710905643e-23, 1.706258568414198e-39, 1.6382855298440747e-20]
    zs_eps_0 = [0.13754371891028325, 0.29845155687154623, 0.2546683930289046, 0.08177453852283137, 0.22756179266643456]
    @test abs(Clapeyron.rr_vle_vapor_fraction(Ks_eps_0,zs_eps_0)) < eps(Float64)
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
    end

    @testset "evalexppoly" begin
        @test Clapeyron.evalexppoly(2,(1,2,3),(3,2,1)) == 22
    end
end
@printline