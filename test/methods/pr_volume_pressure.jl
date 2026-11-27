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

@testset "Constant-shape cubic volume pressure and translation" begin
    binary = ["propane", "butane"]
    ternary = [binary; "pentane"]
    shifts = [3e-6, -8e-6, 1.4e-5] # m³/mol, both signs
    shifted = PR(ternary; translation=ConstantTranslation,
                 translation_userlocations=(; v_shift=shifts))
    plain = PR(ternary)
    peneloux = PR(binary; translation=PenelouxTranslation)
    pr2 = PR(binary)
    for (model, untranslated, T, x) in (
            (shifted, plain, 300., [0.26, 0.31, 0.43]),
            (peneloux, pr2, 280., [0.43, 0.57]),
            (SRK(binary), nothing, 290., [0.43, 0.57]),
            (vdW(binary), nothing, 280., [0.43, 0.57]))
        for side in (Val(:bubble), Val(:dew))
            bubble = side === Val(:bubble)
            pub, plain_method, opt, f!, z, pivot, active = pr_volume_case(model, T, x, side)
            @test VPR._select_volume_pressure_method(opt, model, T, x) === opt
            n = length(x)
            lb = [-Inf; -Inf; zeros(n-1)]; ub = [Inf; Inf; ones(n-1)]
            seed = copy(z)
            @test VPR._pr_volume_try(f!, model, T, x, z, lb, ub, plain_method, active,
                                     pivot, bubble, Val(n)) !== nothing
            @test z == seed
            candidate = pub(model, T, x, opt)
            @test pr_volume_close(candidate, pub(model, T, x, plain_method))
            pr_volume_physical(model, T, x, candidate, bubble)
            if untranslated !== nothing
                c = VPR.translation(model, 1., T, x)
                p, vl, vv, w = candidate
                xl, yv = bubble ? (x,w) : (w,x)
                @test pr_volume_close((p, vl+dot(c,xl), vv+dot(c,yv), w),
                                      pub(untranslated, T, x, plain_method))
                @test abs(dot(c,xl)-dot(c,yv)) > 1e-8
                # The pressure row shear is not zero at a perturbed seed.
                z[1] += 0.06
                physical = f!(similar(z), z)
                shear = c .* (physical[end]*VPR.p_scale(model,x)/(VPR.Rgas(model)*T))
                rectified = [physical[1:n] .+ shear; physical[end]]
                @test norm(shear, Inf) > 1e-5
                obj = VPR._pr_volume_objective(f!, model, T, x, pivot, Val(n), bubble)
                @test pr_volume_err(obj.F(similar(z), z), rectified) < 1e-11
            end
        end
    end
    for (method,wrap) in ((ChemPotBubblePressure(),VPR._PRVolumeBubblePressure),
                          (ChemPotDewPressure(),VPR._PRVolumeDewPressure))
        fast=wrap(method); x=[0.43,0.57]
        @test VPR._select_volume_pressure_method(fast,RKPR(binary),290.,x) === fast
        @test VPR._select_volume_pressure_method(fast,PR(binary;translation=VPR.MTTranslation),290.,x) === method
        @test VPR._select_volume_pressure_method(fast,shifted,290.,[0.43,0.57,0.]) === method
    end
end

@testset "Composition-dependent Delta cubic volume pressure" begin
    names = ["propane", "butane"]
    x = [.45,.55]
    T = 295.
    for eos in (RKPR,Clausius,PatelTeja)
        model = eos(names)
        for side in (Val(:bubble),Val(:dew))
            pr_volume_test(model,T,x,side)
        end
        shifts = [3e-6,-8e-6]
        # A translation-invariance check must keep the fitted EOS fixed.
        # RKPR's default constructor can refit k1 through physical lb_volume.
        alpha_data = if model isa RKPR
            params = model.alpha.params
            (;k1=copy(params.k1.values),k2=copy(params.k2.values),
               acentricfactor=copy(params.acentricfactor.values))
        else
            String[]
        end
        translated = eos(names; alpha_userlocations=alpha_data, translation=ConstantTranslation,
                         translation_userlocations=(;v_shift=shifts))
        @test VPR.α_function(model,1.,T,x,model.alpha) ==
              VPR.α_function(translated,1.,T,x,translated.alpha)
        @test model.params.a.values == translated.params.a.values
        @test model.params.b.values == translated.params.b.values
        for (side,public_call) in ((Val(:bubble),bubble_pressure),(Val(:dew),dew_pressure))
            bubble = side === Val(:bubble)
            pub,stock,opt,f!,seed,pivot,active = pr_volume_case(translated,T,x,side)
            @test VPR._select_volume_pressure_method(opt,translated,T,x) === opt
            lb = [-Inf,-Inf,0.]; ub = [Inf,Inf,1.]
            z0 = copy(seed)
            @test VPR._pr_volume_try(f!,translated,T,x,seed,lb,ub,stock,active,pivot,bubble,Val(2)) !== nothing
            @test seed==z0
            cand = pub(translated,T,x,opt)
            @test pr_volume_close(cand,pub(translated,T,x,stock))
            pr_volume_physical(translated,T,x,cand,bubble)
            p,vl,vv,w=cand
            xl,yv = bubble ? (x,w) : (w,x)
            @test pr_volume_close((p,vl+dot(shifts,xl),vv+dot(shifts,yv),w),public_call(model,T,x,stock))
            @test abs(dot(shifts,xl)-dot(shifts,yv)) > 1e-8
            q=copy(seed);q[1]+=0.06
            function rectified(v)
                physical=f!(similar(v),v)
                shear=shifts .* (physical[end]*VPR.p_scale(translated,x)/(VPR.Rgas(translated)*T))
                return [physical[1:2] .+ shear;physical[end]]
            end
            obj=VPR._pr_volume_objective(f!,translated,T,x,pivot,Val(2),bubble)
            F=similar(q);J=zeros(3,3)
            obj.FJ(F,J,q)
            @test pr_volume_err(F,rectified(q)) < 2e-11
            @test pr_volume_err(J,ForwardDiff.jacobian(rectified,q)) < 2e-10
            @test pr_volume_err(J,pr_volume_fd(rectified,q,pivot)) < 2e-5
        end
    end
end

@testset "Smooth real denominator at coinciding and conjugate Delta roots" begin
    names=["propane","butane"]
    # A genuine RKPR stock model at a synthetic root crossing, not a pole.
    rk=RKPR(names;
        userlocations=(;Tc=[369.83,425.12],Pc=[4.21e6,3.77e6],Vc=[0.000203,0.000255],
            Mw=[44.09562,58.1222],acentricfactor=[0.149,0.197],delta=[0.2,0.7],
            k=zeros(2,2),l=zeros(2,2)),
        alpha_userlocations=(;k1=[0.95,1.08],k2=[2.,2.],acentricfactor=[0.149,0.197]))
    q=[0.00025,295.,0.5287128712871286,0.4712871287128714]
    @test isequal(VPR.cubic_Δ(rk,q[3:4])...)
    f(v)=VPR.a_res(rk,v[1],v[2],v[3:4])
    # Independent 80-decimal integration of the rational kernel and its full derivatives.
    refg=[4851.819265378636,0.007907442627940985,-0.9063944635785931,-1.5568691616746572]
    refH=[-30261118.23351308 -28.068285243117735 1816.0469125910646 3720.2223626176516;
        -28.068285243117735 -5.675927476953798e-5 0.005462504104937247 0.008761060571114995;
        1816.0469125910646 0.005462504104937247 0.3386538840473714 0.5799693039996616;
        3720.2223626176516 0.008761060571114995 0.5799693039996616 0.6793657529298895]
    @test f(q) ≈ -1.4856967973443785 rtol=5e-13
    @test all(isapprox.(ForwardDiff.gradient(f,q),refg;rtol=5e-12,atol=5e-13))
    @test all(isapprox.(ForwardDiff.hessian(f,q),refH;rtol=5e-12,atol=5e-13))
    # A rational-kernel quadrature independent of the log/atan/series formulas.
    # Transform Gauss-Legendre nodes from [-1,1] to [0,1].
    n=32
    ev=eigen(SymTridiagonal(zeros(n),[i/sqrt(4i^2-1) for i in 1:n-1]))
    nodes=(ev.values .+ 1)./2; weights=abs2.(ev.vectors[1,:])
    pt=PatelTeja(names); T=295.; V=.00028; z=[.45,.55]
    c=VPR.diagvalues(pt.params.c.values)
    b=VPR.cubic_lb_volume(pt,T,z)
    γ0=-3+2sqrt(2)
    # Hold composition fixed while moving shape from two real factors,
    # through the repeated limit, to complex conjugate factors.
    for γ in (γ0+0.01,γ0+1e-9,γ0,γ0-1e-9,γ0-0.01,-.4)
        c .= γ*b .+ [.18,-.18*z[1]/z[2]].*b
        _,Δ2=VPR.cubic_ΔT(pt,T,z)
        if γ<γ0-1e-12
            @test !iszero(imag(Δ2))
        elseif γ>γ0+1e-12
            @test iszero(imag(Δ2))
        end
        function rational(v)
            vj,tj,zj1,zj2 = v
            zj = [zj1,zj2]
            aj,bj,_ = VPR.cubic_ab(pt,vj,tj,zj)
            nj = zj1+zj2
            xj = bj*nj/vj
            γj = dot(c,zj)/VPR.cubic_lb_volume(pt,tj,zj)
            # D(t)=1+t+γt(1-t) from the unfactored Patel–Teja form.
            Kj = sum(weights[i]/(1+nodes[i]*xj+γj*nodes[i]*xj*(1-nodes[i]*xj)) for i in eachindex(nodes))
            return -log1p(-xj)-aj*nj/(VPR.Rgas(pt)*tj*vj)*Kj
        end
        actual(v)=VPR.a_res(pt,v[1],v[2],v[3:4])
        val=[V,T,z...]
        @test isapprox(actual(val),rational(val);rtol=1e-12)
        @test pr_volume_err(ForwardDiff.gradient(actual,val),ForwardDiff.gradient(rational,val)) < 5e-10
        @test pr_volume_err(ForwardDiff.hessian(actual,val),ForwardDiff.hessian(rational,val)) < 5e-9
        cs=VPR.SVector{2}(z/V)
        coeff=VPR._pr_volume_coeff(pt,T,z,Val(2))
        g,h,p,r,H=VPR._pr_volume_jet(coeff,cs,Val(true))
        ϕ(v)=sum(v)*rational([1.,T,v...])
        @test pr_volume_err(g,ForwardDiff.gradient(ϕ,cs)) < 5e-12
        @test pr_volume_err(H,ForwardDiff.hessian(ϕ,cs)) < 5e-11
        @test pr_volume_err(h,H*cs) < 5e-12
        @test isapprox(p*VPR.Rgas(pt)*T,pressure(pt,V,T,z);rtol=3e-13)
        @test isfinite(r)
    end
end
