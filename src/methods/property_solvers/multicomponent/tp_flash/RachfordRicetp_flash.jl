"""
    RRTPFlash(;numphases = 2;max_steps = 1e4*(numphases-1),time_limit = Inf)

Method to solve non-reactive multicomponent flash problem by using successive substitution in the Rachford-Rice equation.

User must assume a number of phases, `numphases`. If true number of phases is smaller than numphases, model should predict either (a) identical composition in two or more phases, or (b) one phase with negligible total number of moles. If true number of phases is larger than numphases, a thermodynamically unstable solution will be predicted.

The optimizer will stop at `max_steps` evaluations or at `time_limit` seconds
"""
Base.@kwdef struct RRTPFlash <: TPFlashMethod
    numphases::Int = 2
    max_steps::Int = 1e4*(numphases-1)
    tolerance::Float64 = 1e-10
    K0 = nothing
    time_limit::Float64 = Inf
end

#z is the original feed composition, x is a matrix with molar fractions, n is a matrix with molar amounts

function tp_flash_impl(model::EoSModel, p, T, n, method::RRTPFlash)
    numspecies = length(model)
    
    if method.K0===nothing
        pure = split_model.(model)
        crit = crit_pure.(pure)
        Tc = [crit[i][1] for i ∈ @comps]
        pc = [crit[i][2] for i ∈ @comps]
        ω = acentric_factor.(pure)
        K0 = @. exp(log(pc/p)+5.373*(1+ω)*(1-Tc/T))
    end
    tol = 1.
    
    x = []
    y = []
    α₀ = 0.
    iter = 0
    while tol>method.tolerance && iter<method.max_steps
        g(α) = sum(n[i]*(K0[i]-1.)/(K0[i]*α+(1. -α)) for i ∈ @comps)
        α₀ = Roots.find_zero(g, 0.5)
        
        x = @.  n / (K0*α₀+(1-α₀))
        y = @.  n*K0 / (K0*α₀+(1-α₀))
        φ_α = fugacity_coefficient(model,p,T,x)
        
        φ_β = fugacity_coefficient(model,p,T,y)
        K1 = φ_α./φ_β
        tol = (sum((1-K1[i]/K0[i])^2 for i ∈ @comps))^0.5
        K0 = deepcopy(K1)
        iter += 1
    end

    G = (gibbs_free_energy(model,p,T,x)*(1-α₀)+gibbs_free_energy(model,p,T,y)*α₀)/R̄/T

    X = hcat(x,y)'
    nvals = X.*[1-α₀
                α₀]
    return (X, nvals, G)
end

export RRTPFlash