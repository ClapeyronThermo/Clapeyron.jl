using Test





@testset "RachfordRice.jl — API" begin

# ------------------------------------------------------------------
@testset "RachfordRiceResult{Float64}() — failure placeholder" begin
    r = RachfordRiceResult{Float64}()
    @test isnan(r.V)
    @test isnan(r.L)
end

# ------------------------------------------------------------------
@testset "RachfordRiceResidual — construction & container protocol" begin
    K = [2.0, 1.5, 0.5, 0.1]
    z = [0.25, 0.25, 0.25, 0.25]
    RR = RachfordRiceResidual(K, z)
    @test length(RR) == 4
    @test size(RR)   == (4,)
    @test eltype(RR) == Float64
    @test RR[2]      == 1.5
    @test collect(RR) == K
    @test firstindex(RR) == 1
    @test lastindex(RR)  == 4
    @test eachindex(RR) == 1:4
    @test collect(Iterators.take(RR, 2)) == [2.0, 1.5]

    # length mismatch → BoundsError
    @test_throws BoundsError RachfordRiceResidual([1.0, 2.0], [0.5])
end

@testset "RachfordRiceResidual — non_inx / non_iny masking" begin
    K  = [2.0, 3.0, 0.5]
    z  = [0.3, 0.3, 0.4]
    nx = [true,  false, false]         # component 1 non-condensable  → K → Inf
    ny = [false, true,  false]         # component 2 non-volatile     → K → 0
    RR = RachfordRiceResidual(K, z; non_inx = nx, non_iny = ny)
    @test RR[1] == Inf
    @test RR[2] == 0.0
    @test RR[3] == 0.5
    @test collect(RR) == [Inf, 0.0, 0.5]

    @test_throws BoundsError RachfordRiceResidual(K, z; non_inx = [true])
end

# ------------------------------------------------------------------
@testset "rachfordrice — both signatures give the same answer" begin
    K = [2.0, 1.5, 0.5, 0.1]
    z = [0.25, 0.25, 0.25, 0.25]
    r1 = rachfordrice(K, z)
    r2 = rachfordrice(RachfordRiceResidual(K, z))
    @test r1.V ≈ r2.V atol=1e-14
    @test r1.L ≈ r2.L atol=1e-14
    @test r1.status == r2.status
    @test r1.solved == r2.solved
end

@testset "rachfordrice — user-supplied β0 is honoured when in bounds" begin
    K = [2.0, 1.5, 0.5, 0.1]
    z = [0.25, 0.25, 0.25, 0.25]
    RR = RachfordRiceResidual(K, z)
    lo, hi = rr_limits(RR)
    r = rachfordrice(RR, (lo + hi) / 2)
    @test r.status == RREq
    @test r.solved
end

@testset "rachfordrice — z need not be normalised" begin
    K = [2.0, 1.5, 0.5, 0.1]
    z = [1.0, 1.0, 1.0, 1.0]              # un-normalised
    zn = z ./ sum(z)
    r  = rachfordrice(K, z)
    rn = rachfordrice(K, zn)
    @test r.V ≈ rn.V atol=1e-12
    @test r.status == rn.status
end

# ------------------------------------------------------------------
@testset "rachfordrice — single-phase and trivial results" begin
    # all K < 1 → single-phase liquid
    r = rachfordrice([0.5, 0.3], [0.5, 0.5])
    @test r.status == RRLiquid
    @test r.V == 0.0
    @test r.L == 1.0
    @test r.solved == false
    @test r.feasible == false

    # all K > 1 → single-phase vapour
    r = rachfordrice([2.0, 3.0], [0.5, 0.5])
    @test r.status == RRVapour
    @test r.V == 1.0
    @test r.L == 0.0
    @test r.solved == false
    @test r.feasible == false

    z1 = [0.0, 0.0885053990596404, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.721469918219507, 0.01742948033685831,
            0.1725952023839942]
    K1 = [192.3625321718047, 105.20070698573475, 76.30532397111489, 37.090890262982, 21.38862102676539, 18.093547012968767,
            10.319129068837443, 9.001962137200403, 4.565198737490148, 340.69153314749224, 269.09234343328467,
            21.858385052861507]
    r = rachfordrice(K1,z1)
    @test r.V == 1.0
    @test r.L == 0.0
    @test r.solved == false
    @test r.feasible == false

    # all K == 1 → trivial
    r = rachfordrice([1.0, 1.0], [0.5, 0.5])
    @test r.status == RRTrivial
    @test isnan(r.V)
    @test r.solved == false

    # negative K → failure
    r = rachfordrice([-1.0, 2.0], [0.5, 0.5])
    @test r.status == RRFailure
    @test isnan(r.V)
end

@testset "rachfordrice — negative_flash for a β>1 root" begin
    # constructed so that the two-phase root is at β = 8.9 (liquid-like)
    K = [2.0, 0.9]; z = [0.9, 0.1]
    r_phys   = rachfordrice(K, z)                          # → (1, 0), RRVapour
    r_neg    = rachfordrice(K, z; negative_flash = true)   # → β ≈ 8.9
    @test r_phys.status == RRVapour
    @test r_phys.V == 1.0
    @test r_phys.solved == false
    @test r_neg.status == RRVapour
    @test r_neg.negative_flash == true
    @test r_neg.solved == true
    @test r_neg.V ≈ 8.9 atol=1e-6
    # residual at the negative-flash root is ~0
    RR = RachfordRiceResidual(K, z)
    @test abs(residual(RR, r_neg.V)) ≤ 1e-12
end

# ------------------------------------------------------------------
@testset "residual family" begin
    K = [2.0, 1.5, 0.5, 0.1]; z = fill(0.25, 4)
    RR = RachfordRiceResidual(K, z)
    β  = 0.4
    F = residual(RR, β)
    F2, dF = residual_first_derivative(RR, β)
    F3, dF3, d2F = residual_first_second_derivative(RR, β)
    F4, dF4, d2F4, S = residual_first_second_derivative_error(RR, β)

    @test F == F2 == F3 == F4
    @test dF == dF3 == dF4
    @test d2F == d2F4
    @test S > 0
    # βl default
    @test residual(RR, β, 1 - β) == F
    # sign of derivative
    @test dF ≤ 0
end

# ------------------------------------------------------------------
@testset "rr_limits / rr_status" begin
    K = [2.0, 1.5, 0.5, 0.1]; z = fill(0.25, 4)
    RR = RachfordRiceResidual(K, z)
    lo, hi = rr_limits(RR)
    @test lo < hi
    r = rachfordrice(RR)
    @test lo ≤ r.V ≤ hi
    @test rr_status(RR) == RREq

    # single-phase → NaN limits
    RRl = RachfordRiceResidual([0.4, 0.2], [0.5, 0.5])
    lo, hi = rr_limits(RRl)
    @test isnan(lo) && isnan(hi)
end

# ------------------------------------------------------------------
@testset "vapour_compositions / liquid_compositions — all signatures" begin
    K = [2.0, 1.5, 0.5, 0.1]; z = fill(0.25, 4)
    RR = RachfordRiceResidual(K, z)
    res = rachfordrice(RR)
    βv, βl = res.V, res.L

    y1 = vapour_compositions(K, z, βv)               # default βl
    y2 = vapour_compositions(K, z, βv, βl)
    y3 = vapour_compositions(RR, βv)
    y4 = vapour_compositions(RR, βv, βl)
    y5 = vapour_compositions(RR, res)                # RachfordRiceResult overload
    @test y1 == y2 == y3 == y4 == y5

    x1 = liquid_compositions(K, z, βv)
    x2 = liquid_compositions(RR, βv)
    x3 = liquid_compositions(RR, res)
    @test x1 == x2 == x3

    ybuf = similar(z); xbuf = similar(z)
    vapour_compositions!(ybuf, K, z, βv, βl);      @test ybuf == y1
    vapour_compositions!(ybuf, RR, βv);            @test ybuf == y1
    liquid_compositions!(xbuf, K, z, βv, βl);      @test xbuf == x1
    liquid_compositions!(xbuf, RR, βv);            @test xbuf == x1

    @test sum(y1) ≈ 1 atol=1e-12
    @test sum(x1) ≈ 1 atol=1e-12
    @test βv .* y1 .+ βl .* x1 ≈ z atol=1e-12
end

# ------------------------------------------------------------------
@testset "RachfordRiceConvergence — basic behaviour" begin
    x_ok  = [0.3, 0.3, 0.4]
    x_bad = [0.3, 0.3, 0.5]
    s1, ok1 = Conv.composition_convergence_test(x_ok)
    s2, ok2 = Conv.composition_convergence_test(x_bad)
    @test ok1 && s1 ≤ 0
    @test !ok2 && s2 > 0

    s, ok = Conv.phase_fraction_convergence_test(0.4, 0.6)
    @test ok && s ≤ 0
    s, ok = Conv.phase_fraction_convergence_test(0.4, 0.7)
    @test !ok && s > 0

    # "good" RR solution converges
    K = [2.0, 1.5, 0.5, 0.1]; z = fill(0.25, 4)
    RR = RachfordRiceResidual(K, z)
    r  = rachfordrice(RR)
    y  = vapour_compositions(RR, r)
    x  = liquid_compositions(RR, r)
    all_ok, sev, cv = Conv.check_convergence(y, x, r.V, r.L, z, K)
    @test all_ok
    @test sev.all_converged
    @test cv.within_bounds
end

end # testset "API"