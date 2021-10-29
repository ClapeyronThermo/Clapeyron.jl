function Tp_flash(model::EoSModel,T,p,z;phases=:VLE,K0=nothing)
    nc = length(z)
    if K0 == nothing 
        # Add initial guesses function
    end
    tol = 1.
    # println(K0)
    x = []
    y = []
    α₀ = 0.
    while tol>1e-10
        g(α) = sum(z[i]*(K0[i]-1.)/(K0[i]*α+(1. -α)) for i ∈ @Clapeyron.comps)
        α₀ = find_zero(g, 0.5)
        
        x = @.  z / (K0*α₀+(1-α₀))
        y = @.  z*K0 / (K0*α₀+(1-α₀))
#         println(x)
#         println(y)
        φ_α = fugacity_coefficient(model,p,T,x)
        if phases==:VLE
            φ_β = fugacity_coefficient(model,p,T,y)
        else
            φ_β = fugacity_coefficient(model,p,T,y;phase=:liquid)
        end
        K1 = φ_α./φ_β
        tol = (sum((1-K1[i]/K0[i])^2 for i ∈ @Clapeyron.comps))^0.5
        K0 = deepcopy(K1)
#         println(φ_α)
#         println(φ_β)
    end
    return x,y,α₀
end

function Tp_flash(model::ElectrolyteModel,T,p,z;phases=:VLE,K0=nothing)
    nc = length(z)
    if K0 == nothing 
        # Add initial guesses function
    end
    tol = 1.
    # println(K0)
    x = []
    y = []
    α₀ = 0.
    while tol>1e-3
        g(α) = sum(z[i]*(K0[i]-1.)/(K0[i]*α+(1. -α)) for i ∈ @Clapeyron.comps)
        α₀ = find_zero(g, 0.5)
        
        x = @.  z / (K0*α₀+(1-α₀))
        y = @.  z*K0 / (K0*α₀+(1-α₀))
        φ_α = fugacity_coefficient(model,p,T,x)
        if phases==:VLE
            φ_β = fugacity_coefficient(model,p,T,y)
        else
            φ_β = fugacity_coefficient(model,p,T,y;phase=:liquid)
        end
        K1 = φ_α./φ_β
        ν = model.stoic_coeff
        if phases == :LLE
            K1[end-1:end] .= (K1[end-1]^ν[1][1]*K1[end]^ν[1][2])^(1/(ν[1][1]+ν[1][2]))
        else
            K1[end-1:end] .= 1e-30
        end
        
        tol = (sum((1-K1[i]/K0[i])^2 for i ∈ @Clapeyron.comps))^0.5
        K0 = deepcopy(K1)
#         println(tol)
#         println(φ_α)
#         println(φ_β)
    end

    return x,y,α₀
    end
