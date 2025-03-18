struct MSAIDParam <: EoSParam
    sigma::SingleParam{Float64}
    dipole::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type MSAIDModel <: IonModel end

struct MSAID <: MSAIDModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::MSAIDParam
    references::Array{String,1}
end

@registermodel MSAID

export MSAID
"""
    MSAID(solvents::Array{String,1}, 
         ions::Array{String,1}; 
         SAFTlocations=String[], 
         userlocations=String[], 
         verbose=false)

## Input parameters
- `sigma`: Single Parameter (`Float64`) - Hard-sphere diameter `[m]`
- `charge`: Single Parameter (`Float64`) - Charge `[-]`

## Description
This function is used to create a Mean Spherical Approximation model. The MSAID term gives the excess Helmholtz energy to account for the electrostatic interactions between ions in solution.

## References
1. Blum, L. (1974). Solution of a model for the solvent‐electrolyte interactions in the mean spherical approximation, 61, 2129–2133.
"""
function MSAID(solvents,ions; userlocations, verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)
    icomponents = 1:length(components)
    params = getparams(components, ["Electrolytes/properties/charges.csv","properties/molarmass.csv"]; userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    if any(keys(params).=="b")
        params["b"].values .*= 3/2/N_A/π*1e-3
        params["b"].values .^= 1/3
        sigma = SingleParam("sigma",components,params["b"].values)
    else
        params["sigma"].values .*= 1E-10
        sigma = params["sigma"]
    end

    dipole = params["dipole"]
    dipole.values .*= 1/(299792458)*1e-21
    
    charge = params["charge"]

    packagedparams = MSAIDParam(sigma,dipole,charge)

    references = String[]
        
    model = MSAID(components, icomponents, packagedparams, references)
    return model
end

function a_res(model::MSAIDModel, V, T, z, _data=@f(data))
    return a_ion(model, V, T, z, _data)
end

function data(model::MSAIDModel, V, T, z)
    β = 1/(k_B*T)
    σ = model.params.sigma.values
    Z = model.params.charge.values

    isolv = model.icomponents[Z.==0]
    iions = model.icomponents[Z.!=0]

    μ = model.params.dipole.values[isolv][1]

    ρ = N_A*sum(z)/V
    x = z/sum(z)
    
    σₙ = σ[isolv][1]
    xₙ = x[isolv][1]
    ρₙ = ρ*xₙ

    α₀ = e_c*√(β/ϵ_0) # Checked
    α₂ = μ*√(β/3/ϵ_0) # Checked
    Δ = 1-π*ρ/6*sum(x[i]*σ[i]^3 for i ∈ @comps)
    ξ₂ = ρ*sum(x[i]*σ[i]^2 for i ∈ @comps)
    χ = sum(x[i]*σ[i]*Z[i] for i in iions)

    _data = (α₀,α₂,Δ,ξ₂,χ,σₙ,ρₙ,ρ,x)
    return _data
end

function obj_MSAID(F,model::MSAIDModel,Γ,B,b₂,_data)
    # println(Γ)
    # println(B)
    # println(b₂)
    nc = length(model)
    σ = model.params.sigma.values
    Z = model.params.charge.values

    isolv = model.icomponents[Z.==0]
    iions = model.icomponents[Z.!=0]

    (α₀,α₂,Δ,ξ₂,χ,σₙ,ρₙ,ρ,x) = _data

    β₆ = 1-b₂/6 # Checked
    λ  = (1+b₂/3)/β₆ # Checked
    y₁ = 4/(β₆*(1+λ)^2) # Checked

    σΓ = σ*Γ
    W₁ = ρ*sum(x[i]*Z[i]^2/(β₆*(σₙ+σ[i]*λ)*(1+σΓ[i])) for i in iions) # Checked
    W₂ = 1/2*ρₙ*ρ*σₙ^2*B*sum(x[i]*σ[i]^2*Z[i]^2/(2*β₆*(σₙ+σ[i]*λ)*(1+σΓ[i]))^2 for i in iions) # Checked
    Vη = (-W₁/2+√((W₁/2)^2+2B*W₂/β₆^2))/(W₂) # Checked

    ΔΓ = @. Vη*ρₙ*σₙ^2*σ^2*B/(8*β₆*(σₙ+λ*σ)) # Checked
    Dᶠ = @. Z*β₆/(2*(1+σΓ-ΔΓ)) # Checked

    D = 1 + Vη^2*ρₙ*ρ*σₙ^2*sum(x[i]*σ[i]^2*Dᶠ[i]^2/(2β₆*(σₙ+λ*σ[i]))^2 for i in iions) # Checked

    Dac = ρ*sum(x[i]*Dᶠ[i]^2 for i in iions) # Checked
    Ω   = Vη*ρ*sum(x[i]*σ[i]*Dᶠ[i]^2/(σₙ+λ*σ[i]) for i in iions) # Checked

    Γₛ  = @.  ((1+σΓ-ΔΓ)*D-1)/σ # Checked

    a⁰  = @. β₆*Γₛ*Dᶠ/Dac # Checked
    a¹  = D*β₆*(σₙ*B/2 + Ω*λ/(D*β₆))/(2*Dac) # Checked
    
    K¹⁰ = @. -(σₙ^2*Dᶠ*(Vη/(σₙ+λ*σ)+Ω*Γₛ/Dac)/(2*D*β₆^2) + σₙ^3*B*a⁰/(12*β₆)) # Checked
    K¹¹ = (1-(λ+ρₙ*σₙ^2*Ω*a¹/(2*β₆^2))/(D*β₆)-ρₙ*σₙ^3*B*a¹/(12*β₆))/ρₙ # Checked
    
    F[1] = (ρ*sum(x[i]*a⁰[i]^2 for i in iions) + ρₙ*a¹^2)/α₀^2 - 1 # Checked
    F[2] = (-ρ*sum(x[i]*a⁰[i]*K¹⁰[i] for i in iions) + a¹*(1-ρₙ*K¹¹))/(α₀*α₂) - 1 # Checked
    F[3] = ((1-ρₙ*K¹¹)^2+ρₙ*ρ*sum(x[i]*K¹⁰[i]^2 for i in iions) - y₁^2)/(ρₙ*α₂^2) - 1 # Checked

    η = ρ*sum(Z.^2 .*x)
    m = @. Vη*Dᶠ/(σₙ+λ*σ) * √(η*ρₙ)*σₙ*σ/Z
    N = @. (2*Dᶠ/(β₆*σ)*(1+Vη*ρₙ*σₙ^3*B*σ/(24*(σₙ+λ*σ))) - Z/σ)*σ/Z
    ϵr = 1+ρₙ*α₂^2*β₆^2*(1+λ)^4/16
    return F, m, N, ϵr
end

function a_ion(model::MSAIDModel, V, T, z, _data=@f(data))
    σ = model.params.sigma.values
    Z = model.params.charge.values

    isolv = model.icomponents[Z.==0]
    iions = model.icomponents[Z.!=0]

    nc = length(model)
    

    isolv = model.icomponents[Z.==0]
    iions = model.icomponents[Z.!=0]

    (α₀,α₂,Δ,ξ₂,χ,σₙ,ρₙ,ρ,x) = _data

    Γ, B, b₂ = @f(solve_MSAID, _data)

    β₆ = 1-b₂/6 # Checked
    λ  = (1+b₂/3)/β₆ # Checked
    y₁ = 4/(β₆*(1+λ)^2) # Checked

    σΓ = σ*Γ
    W₁ = ρ*sum(x[i]*Z[i]^2/(β₆*(σₙ+σ[i]*λ)*(1+σΓ[i])) for i in iions) # Checked
    W₂ = 1/2*ρₙ*ρ*σₙ^2*B*sum(x[i]*σ[i]^2*Z[i]^2/(2*β₆*(σₙ+σ[i]*λ)*(1+σΓ[i]))^2 for i in iions) # Checked
    Vη = (-W₁/2+√((W₁/2)^2+2B*W₂/β₆^2))/(W₂) # Checked

    ΔΓ = @. Vη*ρₙ*σₙ^2*σ^2*B/(8*β₆*(σₙ+λ*σ)) # Checked
    Dᶠ = @. Z*β₆/(2*(1+σΓ-ΔΓ)) # Checked

    D = 1 + Vη^2*ρₙ*ρ*σₙ^2*sum(x[i]*σ[i]^2*Dᶠ[i]^2/(2β₆*(σₙ+λ*σ[i]))^2 for i in iions) # Checked

    Dac = ρ*sum(x[i]*Dᶠ[i]^2 for i in iions) # Checked
    Ω   = Vη*ρ*sum(x[i]*σ[i]*Dᶠ[i]^2/(σₙ+λ*σ[i]) for i in iions) # Checked

    Γₛ  = @.  ((1+σΓ-ΔΓ)*D-1)/σ # Checked

    a¹  = D*β₆*(σₙ*B/2 + Ω*λ/(D*β₆))/(2*Dac) # Checked

    Jp = zero(one(eltype(model))+V+T+first(z))

    for i in iions
        σi = σ[i]
        for j in iions
            σij = (σi+σ[j])/2
            Q00 = 2π/Δ*(σij+π*σi*σ[j]*ξ₂/(4Δ)) 
                - 1/2*Dᶠ[i]*Dᶠ[j]*(ρₙ*σₙ^2*Vη^2/(D*β₆^2*(σₙ+λ*σi)*(σₙ+λ*σ[j]))
                + 4*Γₛ[i]*Γₛ[j]/(D*Dac))
            Jp += x[i]*x[j]*σij*Q00^2/(3π)^2
        end
    end

    for i in iions
        σi = σ[i]
        for j in isolv
            σij = (σi+σ[j])/2
            Q00 = 2π/Δ*(σij+π*σi*σ[j]*ξ₂/(4Δ)) 
            Jp += 2*x[i]*x[j]*σij*Q00^2/(3π)^2

            Q01 = Dᶠ[i]/(D*β₆)*(λ*Vη/(σₙ+λ*σi)+2*Γₛ[i]*a¹)
            Jp += 2/9*x[i]*x[j]*σij*Q01^2/(3π)^2
        end
    end

    for i in isolv
        σi = σ[i]
        for j in isolv
            σij = (σi+σ[j])/2
            Q11 = 2λ/(D*ρₙ*σₙ^2)*(λ+ρₙ*σₙ^2*Ω*a¹/(2*β₆))
                + σₙ*B*a¹/(2*β₆) - 2/(ρₙ*σₙ^2)
            
            qp = -b₂*(λ+3)/(1+λ)^2/(ρₙ*σₙ^2)

            h11 = (Q11+2qp)/(2*√(3)*π)
            h12 = √(10)*(Q11-qp)/(2*√(3)*π)
            Jp += x[i]*x[j]*σij*(h11^2 + 1/5*h12^2)
        end
    end

    Jp *= ρ/(3π)

    N = @. (2*Dᶠ/(β₆*σ)*(1+Vη*ρₙ*σₙ^3*B*σ/(24*(σₙ+λ*σ))) - Z/σ)

    return (α₀^2*sum(x[i]*Z[i]*N[i] for i in iions)
    - ρₙ/ρ*α₀*B*α₂)/(6π) - Jp
end

function solve_MSAID(model::MSAIDModel,V,T,z,_data = @f(data))
    F = zeros(typeof(one(eltype(model))+V+T+first(z)),(3))

    x0 = ones(typeof(one(eltype(model))+V+T+first(z)),(3)).*[9., 1.,2.02]

    f!(F,x) = obj_MSAID(F,model,x[1]*1e9,x[2]*1e17,x[3],_data)[1]

    sol = Solvers.nlsolve(f!,x0, LineSearch(Newton()))
    _x = Solvers.x_sol(sol)
    return _x[1]*1e9, _x[2]*1e17, _x[3]
end