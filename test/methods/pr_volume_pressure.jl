using Clapeyron, ForwardDiff, LinearAlgebra, NLSolvers, Test
const VPR = Clapeyron

function pr_volume_case(model, T, x, ::Val{S}) where S
    bubble = S === :bubble
    stock = bubble ? ChemPotBubblePressure() : ChemPotDewPressure()
    opt = bubble ? VPR._PRVolumeBubblePressure(stock) : VPR._PRVolumeDewPressure(stock)
    public_call = bubble ? bubble_pressure : dew_pressure
    initialize = bubble ? VPR.bubble_pressure_init : VPR.dew_pressure_init
    active = trues(length(x))
    _, vl, vv, w = initialize(model, T, x, nothing, nothing, nothing, active, false)
    pivot = argmax(w)
    z = [VPR.η_from_v(model, vl, T, bubble ? x : w),
         VPR.η_from_v(model, vv, T, bubble ? w : x),
         [w[i] for i in eachindex(w) if i != pivot]...]
    function f!(F, q)
        free = @view(q[3:end])
        if bubble
            return VPR.Obj_bubble_pressure(model, model, F, T, q[1], q[2], x, free, active, pivot)
        end
        return VPR.Obj_dew_pressure(model, model, F, T, q[1], q[2], free, x, active, pivot)
    end
    return public_call, stock, opt, f!, z, pivot, active
end

function pr_volume_fd(f, z, pivot)
    J = zeros(length(z), length(z))
    w = VPR.FractionVector(@view(z[3:end]), pivot)
    for j in eachindex(z)
        h = j <= 2 ? 1e-4 : min(1e-4, 0.04 * min(z[j], w[pivot]))
        d = zeros(length(z)); d[j] = h
        J[:, j] .= (f(z - 2d) - 8f(z - d) + 8f(z + d) - f(z + 2d)) / (12h)
    end
    return J
end

pr_volume_err(a, b) = maximum(abs.(a .- b) ./ (1 .+ abs.(b)))

function pr_volume_close(a, b)
    return all(i -> isapprox(a[i], b[i]; rtol=5e-8, atol=2e-11), 1:4)
end

function pr_volume_physical(model, T, x, pvvw, bubble)
    p, vl, vv, w = pvvw
    xl, yv = bubble ? (x, w) : (w, x)
    pl = pressure(model, vl, T, xl); pv = pressure(model, vv, T, yv)
    @test 0 < vl < vv && min(pl, pv, p) > 0
    @test abs(sum(w) - 1) < 1e-13 && all(>(0), w)
    @test abs(pl - pv) / max(pl, pv) < 1e-6
    @test isapprox(p, bubble ? pl : pv; rtol=1e-12)
    lg(v, z, pp) = log.(VPR.VT_fugacity_coefficient(model, v, T, z) .* z .* pp)
    @test norm(lg(vl, xl, pl) - lg(vv, yv, pv), Inf) < 3e-8
    @test VPR.VT_mechanical_stability(model, vl, T, xl) &&
        VPR.VT_mechanical_stability(model, vv, T, yv)
end

function pr_volume_test(model, T, x, side::Val{S}; replay_seeds=false) where S
    pub, stock, opt, f!, z, pivot, active = pr_volume_case(model, T, x, side)
    @test VPR._select_volume_pressure_method(opt, model, T, x) === opt
    bubble = S === :bubble
    lb = [-Inf; -Inf; zeros(length(z)-2)]
    ub = [Inf; Inf; ones(length(z)-2)]
    q0 = copy(z)
    f(q) = f!(similar(q), q)
    F = similar(z); J = zeros(length(z), length(z))
    obj = VPR._pr_volume_objective(f!, model, T, x, pivot, Val(length(x)), bubble)
    for shift in ((0., 0.), (0.04, -0.05), (-0.06, 0.07), (0.02, -1.2))
        q = copy(z)
        q[1] += shift[1]; q[2] += shift[2]
        q[3:end] .+= shift[1] * 0.001 * minimum(VPR.FractionVector(@view(z[3:end]), pivot))
        Ftrial = copy(obj.F(F, q))
        @test pr_volume_err(Ftrial, f(q)) < 2e-11
        obj.FJ(F, J, q)
        @test pr_volume_err(Ftrial, F) < 128eps(Float64)
        @test pr_volume_err(F, f(q)) < 2e-11
        @test pr_volume_err(J, ForwardDiff.jacobian(f, q)) < 2e-10
        @test pr_volume_err(J, pr_volume_fd(f, q, pivot)) < 2e-5
    end
    @test VPR._pr_volume_try(f!, model, T, x, z, lb, ub, stock, active, pivot,
                             bubble, Val(length(x))) !== nothing
    @test z == q0
    baseline = pub(model, T, x, stock)
    @test isequal(pub(model, T, x), baseline)
    candidate = pub(model, T, x, opt)
    @test pr_volume_close(baseline, candidate)
    pr_volume_physical(model, T, x, candidate, bubble)
    if replay_seeds
        for shift in ((0.0, 0.0), (0.02, -0.05), (-0.05, 0.04))
            other = copy(z); other[1] += shift[1]; other[2] += shift[2]
            seed = copy(other)
            reference = VPR.Solvers.nlsolve(f!, copy(other),
                VPR.LineSearch(VPR.Newton2(other), VPR.BoundedLineSearch(lb, ub)),
                NLSolvers.NEqOptions(stock), ForwardDiff.Chunk(min(length(other), 8)))
            result = VPR._pr_volume_solve(f!, model, T, x, other, lb, ub,
                                          stock, pivot, bubble, Val(length(x)))
            @test VPR.__check_convergence(reference) && VPR.__check_convergence(result)
            @test result.options.maxiter == 12
            @test reference.info.iter == result.info.iter
            @test pr_volume_err(reference.info.solution, result.info.solution) < 1e-8
            @test other == seed
        end
    end
end

function pr_volume_route_tests(model, side::Val{S}) where S
    bubble = S === :bubble
    pub = bubble ? bubble_pressure : dew_pressure
    base = bubble ? ChemPotBubblePressure : ChemPotDewPressure
    wrap = bubble ? VPR._PRVolumeBubblePressure : VPR._PRVolumeDewPressure
    m = base()
    x = [0.32, 0.28, 0.40]
    opt = wrap(m)
    @test VPR._select_volume_pressure_method(opt, model, ForwardDiff.Dual(295., 1.), x) === m
    @test VPR._select_volume_pressure_method(opt, model, Float32(295.), x) === m
    @test VPR._select_volume_pressure_method(opt, model, 295., Float32.(x)) === m
    @test VPR._select_volume_pressure_method(opt, model, 295., [1e-30, .4, .6]) === opt
    for z in ([1e-30, .4, .6], [1e-7, .4, .5999999])
        @test isequal(pub(model, 295., z, opt), pub(model, 295., z, m)) ||
            pr_volume_close(pub(model, 295., z, opt), pub(model, 295., z, m))
    end
    dstock = ForwardDiff.derivative(T -> pub(model, T, x, m)[1], 295.)
    dopt = ForwardDiff.derivative(T -> pub(model, T, x, opt)[1], 295.)
    @test isfinite(dstock) && isequal(dstock, dopt)
    @test VPR._select_volume_pressure_method(opt, PCSAFT(["propane", "butane", "pentane"]), 295., x) === m
    @test isequal(pub(model, 295., [0.3, 0.0, 0.7], opt), pub(model, 295., [0.3, 0.0, 0.7], m))
    @test isequal(pub(model, 295., [0.0, 1.0, 0.0], opt), pub(model, 295., [0.0, 1.0, 0.0], m))
    constrained = bubble ? base(nonvolatiles=["pentane"]) : base(noncondensables=["propane"])
    for bypass in (base(max_iters=0), base(atol=1e-9), constrained)
        @test VPR._select_volume_pressure_method(wrap(bypass), model, 295., x) === bypass
        @test isequal(pub(model, 295., x, wrap(bypass)), pub(model, 295., x, bypass))
    end
    for T in (350., 365., 390., 420.)
        stock, candidate = pub(model, T, x, m), pub(model, T, x, opt)
        if isfinite(stock[1])
            @test pr_volume_close(stock, candidate)
        else
            @test isequal(stock, candidate)
        end
    end

    # Force analytic candidate rejection: its true EOS root is the initial
    # state, but the supplied original residual has a different exact zero.
    # This checks both mutation safety and execution of the stock retry, not
    # merely skipping the candidate at the precision/iteration gate.
    p, vl, vv, w = pub(model, 295., x, m)
    @test isfinite(p)
    i = argmax(w)
    q = [VPR.η_from_v(model, vl, 295., bubble ? x : w),
         VPR.η_from_v(model, vv, 295., bubble ? w : x),
         [w[j] for j in eachindex(w) if j != i]...]
    snapshot = copy(q)
    target = copy(q); target[1] += 0.04
    inconsistent!(F, z) = (F .= z .- target)
    lb = [-Inf; -Inf; zeros(length(q)-2)]; ub = [Inf; Inf; ones(length(q)-2)]
    @test VPR._pr_volume_try(inconsistent!, model, 295., x, q, lb, ub, m, trues(3), i,
                             bubble, Val(3)) === nothing
    @test q == snapshot
    candidate, ok = VPR._volume_pressure_solve(VPR._PRVolumeJacobian(), inconsistent!,
        q, lb, ub, m, model, 295., x, trues(3), i, side)
    @test ok && isapprox(candidate, target; rtol=1e-10)
end

@testset "Private exact PR fixed-T volume Jacobian" begin
    pr2 = PR(["propane", "butane"])
    pr3 = PR(["propane", "butane", "pentane"])
    asymmetric = deepcopy(pr3)
    asymmetric.params.b.values[1, 2] *= 1.1
    for side in (Val(:bubble), Val(:dew))
        pr_volume_test(pr2, 295., [0.4, 0.6], side; replay_seeds=true)
        pr_volume_test(pr3, 295., [0.3, 0.3, 0.4], side; replay_seeds=true)
        pr_volume_route_tests(pr3, side)
        @test !VPR._pr_volume_supported(ChemPotBubblePressure(), asymmetric, 295., [0.3, 0.3, 0.4])
    end
    # Nonadditive covolume and nonzero binary interactions, on different
    # species order/composition/temperature than the primary fixtures.
    for (i,j,kij,lij) in ((2,1,.03,.017),(3,2,-.012,.021),(3,1,.019,-.013))
        a=pr3.params.a.values; b=pr3.params.b.values
        a[i,j]=a[j,i]=sqrt(a[i,i]*a[j,j])*(1-kij)
        b[i,j]=b[j,i]=(b[i,i]+b[j,j])*(1-lij)/2
    end
    pr_volume_test(pr3, 305., [0.13, 0.26, 0.61], Val(:bubble))
    pr_volume_test(pr3, 305., [0.13, 0.26, 0.61], Val(:dew))
end
