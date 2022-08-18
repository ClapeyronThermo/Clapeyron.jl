import NLsolve, Optim

function Michelsen_multiphase_tp_flash(model, p, T, z; abstol=1e-10, outerSSiters=30, outerNewtoniters=30)
    w_array, tm_array = chemical_stability_analysis(model, p, T, z; converge_min=true, abstol=1e-3)
    # w_array, tpd_array, _, _ = all_tpd(model, p, T, z)
    # w_array, tpd_array = tpd(model, p, T, z)

    if any(tm_array .< 0.0)
        stable = false

        w_array = w_array[tm_array.<0.0]

        nC = length(z)
        nF = length(w_array)

        # Create initial guesses
        if length(w_array) ≥ 2
            x = zeros(nF, nC)
            β = ones(nF) ./ nF
            for i in 1:nF
                x[i, :] .= w_array[i]
            end
        else # Initialise with K-factors?
            Kʷ = wilson_k_values(model, p, T)
            x[1, :] .= Kʷ
            x[2, :] .= 1 ./ Kʷ
        end
    else
        stable = true
        x = z
        β = [1.0]
    end

    #TODO: Add stability check on resulting phase distribution
    x, β = multiphase_flash_impl(model, p, T, z, x, β, abstol, outerSSiters, outerNewtoniters)

    return x, β
end

function multiphase_E(β, φmat)
    nF, nC = size(φmat)
    E = zeros(eltype(first(β) + first(φmat)), nC)
    for i in range(1, nC)
        for k in range(1, nF)
            E[i] += β[k] / φmat[k, i]
        end
    end
    return E
end

function multiphase_RR_obj(β, z, φmat)
    return sum(β) - sum(z .* log.(multiphase_E(β, φmat)))
end

function multiphase_RR!(β, x, z, φmat)
    nF = length(β)

    fun(x) = multiphase_RR_obj(x, z, φmat)
    function fun_grad!(g, x)
        g .= ForwardDiff.gradient(fun, x)
    end

    function fun_hess!(h, x)
        h .= ForwardDiff.hessian(fun, x)
    end

    # x0 = β
    df = Optim.TwiceDifferentiable(fun, fun_grad!, fun_hess!, β)

    lx = zeros(nF)
    ux = ones(nF)
    dfc = Optim.TwiceDifferentiableConstraints(lx, ux)

    res = Optim.optimize(df, dfc, β, Optim.IPNewton())
    #=
    res = Solvers.optimize(fun,β,NLSolvers.ActiveBox(),bounds = (lx,ux))
    =#
    β .= res.minimizer

    E = multiphase_E(β, φmat)
    # x = zeros(eltype(φmat), size(φmat))
    for (i, xi) in enumerate(eachrow(x))
        xi .= z ./ E ./ φmat[i, :]
    end
end

function fixpoint_multiphase!(F, x, β, φmat, model, p, T, z)
    xold = copy(x)
    multiphase_RR!(β, x, z, φmat)
    # Recalculate fugacity coefficients
    φmat = zeros(eltype(x), size(φmat))
    for (i, φ_i) in enumerate(eachrow(φmat))
        φ_i .= fugacity_coefficient(model, p, T, x[i, :])
    end
    F .= x .- xold
end
NEqOptions(f_limit = method.f_limit,
f_abstol = method.atol,
f_reltol = method.rtol,
maxiter = method.max_iters)
# φ will be an nxm matrix where n (rows) are phases, m (columns) are components
function multiphase_flash_impl(model, p, T, z, x, β, abstol, outerSSiters, outerNewtoniters)
    φmat = zeros(size(x))
    for (i, φ_p) in enumerate(eachrow(φmat))
        φ_p .= fugacity_coefficient(model, p, T, x[i, :])
    end

    f!(F, x) = fixpoint_multiphase!(F, x, β, φmat, model, p, T, z)
    #res = Solvers.nlsolve(f!, x, NLSolvers.Anderson(m=5),NEqOptions(f_abstol = abstol,maxiter = outerSSiters))
    res = NLsolve.nlsolve(f!, x, method=:anderson, m=5, ftol=abstol, iterations=outerSSiters)
    
    if ~(res.x_converged || res.f_converged) # second order scheme
        # @info "SS did not converged in $outerSSiters successive substitution iterations, moving on to second order minimisation"
        #res = Solvers.nlsolve(f!,Solvers.x_sol(res),TrustRegion(Newton(), DogLeg()),NEqOptions(f_abstol = abstol,maxiter = outerNewtoniters))
        res = NLsolve.nlsolve(f!, res.zero, method=:trust_region, autodiff=:forward, ftol=abstol, iterations=outerNewtoniters)
    end
    ~(res.x_converged || res.f_converged) && @warn "Flash calculation failed to converge in $outerSSiters successive substitution iterations and $outerNewtoniters Newton iterations"

    return res.zero, β
end

