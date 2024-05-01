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

    lx = zeros(nF)
    ux = ones(nF)

    fact(x) = cholesky(Solvers.PositiveFactorizations.Positive,x) 
    res1 = Solvers.optimize(fun,β,NLSolvers.ActiveBox(factorize = fact),bounds = (lx,ux))
    β .= Solvers.x_sol(res1)

    E = multiphase_E(β, φmat)
    # x = zeros(eltype(φmat), size(φmat))
    for (i, xi) in enumerate(eachrow(x))
        xi .= z ./ E ./ φmat[i, :]
    end
end

function fixpoint_multiphase!(F, x, β, φmat, model, p, T, z)
    sφ = size(φmat) #at the moment, NLSolvers doesn't support arbitrary shapes :( https://github.com/JuliaNLSolvers/NLSolvers.jl/issues/25
    x = reshape(x,sφ)
    F = reshape(F,sφ)
    xold = copy(x)
    multiphase_RR!(β, x, z, φmat)
    # Recalculate fugacity coefficients
    φmat = zeros(eltype(x), size(φmat))
    for (i, φ_i) in enumerate(eachrow(φmat))
        φ_i .= fugacity_coefficient(model, p, T, x[i, :])
    end
    F .= x .- xold
    return vec(F)
end

# φ will be an nxm matrix where n (rows) are phases, m (columns) are components
function multiphase_flash_impl(model, p, T, z, x, β, abstol, outerSSiters, outerNewtoniters)
    φmat = zeros(size(x))
    for (i, φ_p) in enumerate(eachrow(φmat))
        φ_p .= fugacity_coefficient(model, p, T, x[i, :])
    end

    f!(F, x) = fixpoint_multiphase!(F, x, β, φmat, model, p, T, z)
    method = AndersonFixPoint(;picard_damping=0.5,damping=0.5,memory=5,delay=5,drop_tol=Inf)

    res = Solvers.fixpoint(f!,x, AndersonFixPoint(),maxiters = outerSSiters,atol = abstol,return_last = true)
    
    f!(x,res)

    

    #neq_opts = NEqOptions(f_abstol = abstol,maxiter = outerSSiters)
    #res1 = Solvers.nlsolve(f!, x, NLSolvers.Anderson(0,5,1.0,1e-10),neq_opts)
    #res = NLsolve.nlsolve(f!, x, method=:anderson, m=5, ftol=abstol, iterations=outerSSiters)
    #@show res1
    #@show res
    #@show Solvers.x_sol(res1) - vec(res.zero)
    #@show Solvers.x_sol(res1) - vec(x)
    #@show res.f_converged
    #@show res1
    #if !Solvers.converged(res)
    if Solvers.rtol_anderson(x,res) >= abstol # second order scheme
        # @info "SS did not converged in $outerSSiters successive substitution iterations, moving on to second order minimisation"
        newton_res = Solvers.nlsolve(f!, x, TrustRegion(Newton(), Dogleg()), NEqOptions(f_abstol = abstol,maxiter = outerNewtoniters))
        
        #res1 = Solvers.nlsolve(f!,vec(res.zero),TrustRegion(Newton(), Dogleg()),NEqOptions(f_abstol = abstol,maxiter = outerNewtoniters))
        #res = NLsolve.nlsolve(f!, res.zero, method=:trust_region, autodiff=:forward, ftol=abstol, iterations=outerNewtoniters)
        #@show Solvers.x_sol(res1) - vec(res.zero)
        #@show res
    end
    #~(res.x_converged || res.f_converged) && @warn "Flash calculation failed to converge in $outerSSiters successive substitution iterations and $outerNewtoniters Newton iterations"
    #!Solvers.converged(res) && @warn "Flash calculation failed to converge in $outerSSiters successive substitution iterations and $outerNewtoniters Newton iterations"

    return res.zero, β
end

struct MichelsenMultiTPFlash <: TPFlashMethod end

function tp_flash_impl(model,p,T,z,method::MichelsenMultiTPFlash)
    z, b = Michelsen_multiphase_tp_flash(model, p, T, z; abstol=1e-10, outerSSiters=30, outerNewtoniters=30)
end