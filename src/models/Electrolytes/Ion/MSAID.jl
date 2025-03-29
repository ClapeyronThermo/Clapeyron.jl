struct MSAIDParam <: EoSParam
    sigma::SingleParam{Float64}
    dipole::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type MSAIDModel <: IonModel end

struct MSAID <: MSAIDModel
    components::Array{String,1}
    params::MSAIDParam
    references::Array{String,1}
end

export MSAID
"""
    MSAID(solvents::Array{String,1},
        ions::Array{String,1};
        RSPmodel = nothing,
        userlocations = String[],
        RSPmodel_userlocations = nothing,
        verbose = false)

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

    references = String["10.1063/1.1682224"]
    if count(iszero,charge.values) != 1
        throw(error("MSAID only supports one neutral solvent."))
    end
    model = MSAID(components, packagedparams, references)
    return model
end

IonDependency(ionmodel::MSAIDModel) = IndependentIonModel()

function a_res(model::MSAIDModel, V, T, z, _data=@f(data))
    return a_ion(model, V, T, z, _data)
end

function data(model::MSAIDModel, V, T, z)
    β = 1/(k_B*T)
    σ = model.params.sigma.values
    Z = model.params.charge.values

    nc = length(model)
    isolv = findfirst(iszero,Z)
    μ = model.params.dipole.values[isolv]
    ∑z = sum(z)
    ρ = N_A*∑z/V
    σₙ = σ[isolv]
    xₙ = z[isolv]/∑z
    ρₙ = ρ*xₙ
    α₀ = e_c*√(β/ϵ_0) # Checked
    α₂ = μ*√(β/3/ϵ_0) # Checked
    Δ = 1 - π*ρ/6*sum(z[i]*σ[i]^3 for i ∈ @comps)/∑z
    ξ₂ = ρ*sum(z[i]*σ[i]^2 for i ∈ @comps)/∑z
    χ = sum(z[i]*σ[i]*Z[i] for i ∈ @iions)/∑z
    _data = (α₀,α₂,Δ,ξ₂,χ,σₙ,ρₙ,ρ,∑z)
    return _data
end

function obj_MSAID(model::MSAIDModel,z,Γ,B,b₂,_data)
    # println(Γ)
    # println(B)
    # println(b₂)
    nc = length(model)
    σ = model.params.sigma.values
    Z = model.params.charge.values
    iions = @iions
    (α₀,α₂,Δ,ξ₂,χ,σₙ,ρₙ,ρ,∑z) = _data
    β₆ = 1-b₂/6 # Checked
    λ  = (1+b₂/3)/β₆ # Checked
    y₁ = 4/(β₆*(1+λ)^2) # Checked

    W₁ = ρ*sum(z[i]*Z[i]^2/(β₆*(σₙ + σ[i]*λ)*(1+σ[i]*Γ)) for i ∈ iions)/∑z # Checked
    W₂ = 1/2*ρₙ*ρ*σₙ^2*B*sum(z[i]*σ[i]^2*Z[i]^2/(2*β₆*(σₙ+σ[i]*λ)*(1+σ[i]*Γ))^2 for i ∈ iions)/∑z # Checked
    Vη = (-W₁/2+√((W₁/2)^2+2B*W₂/β₆^2))/(W₂) # Checked


    #[ΔΓ] - > [Dᶠ] -> (D,Dac,Ω) ->  ([Γₛ],[a⁰],a¹)

    D = oneunit(Γ + B + b₂)
    Dac,Ω = zero(D),zero(D)
    D_prefactor = Vη^2*ρₙ*ρ*σₙ^2
    for i ∈ iions
        σᵢ,Zᵢ,zᵢ = σ[i],Z[i],z[i]
        xᵢ = zᵢ/∑z
        ΔΓᵢ = Vη*ρₙ*σₙ*σₙ*σᵢ*σᵢ*B/(8*β₆*(σₙ + λ*σᵢ))
        Dᶠᵢ = Zᵢ*β₆/(2*(1 + σᵢ*Γ - ΔΓᵢ))
        D += D_prefactor*xᵢ*σᵢ*σᵢ*Dᶠᵢ*Dᶠᵢ/(2β₆*(σₙ + λ*σᵢ))^2
        Dac += xᵢ*Dᶠᵢ*Dᶠᵢ
        Ω += xᵢ*σᵢ*Dᶠᵢ*Dᶠᵢ/(σₙ + λ*σᵢ)
    end


    Ω *= Vη*ρ
    Dac *= ρ
    a¹  = D*β₆*(σₙ*B/2 + Ω*λ/(D*β₆))/(2*Dac) # Checked
    K¹¹ = (1-(λ+ρₙ*σₙ^2*Ω*a¹/(2*β₆^2))/(D*β₆)-ρₙ*σₙ^3*B*a¹/(12*β₆))/ρₙ # Checked

    #ΔΓ = @. Vη*ρₙ*σₙ^2*σ^2*B/(8*β₆*(σₙ+λ*σ)) # Checked
    #Dᶠ = @. Z*β₆/(2*(1+σ*Γ-ΔΓ)) # Checked

    #D = 1 + Vη^2*ρₙ*ρ*σₙ^2*sum(x[i]*σ[i]^2*Dᶠ[i]^2/(2β₆*(σₙ+λ*σ[i]))^2 for i in iions) # Checked
    #Dac = ρ*sum(x[i]*Dᶠ[i]^2 for i in iions) # Checked
    #Ω   = Vη*ρ*sum(x[i]*σ[i]*Dᶠ[i]^2/(σₙ+λ*σ[i]) for i in iions) # Checked

    ∑1,∑2,∑3 = zero(D),zero(D),zero(D)
    for i ∈ iions
        σᵢ,Zᵢ,zᵢ = σ[i],Z[i],z[i]
        xᵢ = zᵢ/∑z
        ΔΓᵢ = Vη*ρₙ*σₙ*σₙ*σᵢ*σᵢ*B/(8*β₆*(σₙ + λ*σᵢ))
        Dᶠᵢ = Zᵢ*β₆/(2*(1 + σᵢ*Γ - ΔΓᵢ))
        Γₛᵢ = ((1 + σᵢ*Γ - ΔΓᵢ)*D - 1)/σᵢ
        a⁰ᵢ = β₆*Γₛᵢ*Dᶠᵢ/Dac
        K¹⁰ᵢ = -(σₙ*σₙ*Dᶠᵢ*(Vη/(σₙ+λ*σᵢ)+Ω*Γₛᵢ/Dac)/(2*D*β₆^2) + σₙ^3*B*a⁰ᵢ/(12*β₆))
        ∑1 += xᵢ*a⁰ᵢ*a⁰ᵢ
        ∑2 += xᵢ*a⁰ᵢ*K¹⁰ᵢ
        ∑3 += xᵢ*K¹⁰ᵢ*K¹⁰ᵢ
    end

    F1 = (ρ*∑1 + ρₙ*a¹^2)/α₀^2 - 1
    F2 = (-ρ*∑2 + a¹*(1-ρₙ*K¹¹))/(α₀*α₂) - 1
    F3 = ((1-ρₙ*K¹¹)^2+ρₙ*ρ*∑3 - y₁^2)/(ρₙ*α₂^2) - 1

    return SVector((F1,F2,F3))
    #Γₛ  = @. ((1+σ*Γ-ΔΓ)*D-1)/σ # Checked
    #a⁰  = @. β₆*Γₛ*Dᶠ/Dac # Checked
    #K¹⁰ = @. -(σₙ^2*Dᶠ*(Vη/(σₙ+λ*σ)+Ω*Γₛ/Dac)/(2*D*β₆^2) + σₙ^3*B*a⁰/(12*β₆)) # Checked

    #F[1] = (ρ*sum(x[i]*a⁰[i]^2 for i in iions) + ρₙ*a¹^2)/α₀^2 - 1 # Checked
    #F[2] = (-ρ*sum(x[i]*a⁰[i]*K¹⁰[i] for i in iions) + a¹*(1-ρₙ*K¹¹))/(α₀*α₂) - 1 # Checked
    #F[3] = ((1-ρₙ*K¹¹)^2+ρₙ*ρ*sum(x[i]*K¹⁰[i]^2 for i in iions) - y₁^2)/(ρₙ*α₂^2) - 1 # Checked

    #η = ρ*@sum(Z[i]^2 * x[i])
    #m = @. Vη*Dᶠ/(σₙ+λ*σ) * √(η*ρₙ)*σₙ*σ/Z
    #N = @. (2*Dᶠ/(β₆*σ)*(1+Vη*ρₙ*σₙ^3*B*σ/(24*(σₙ+λ*σ))) - Z/σ)*σ/Z
    #ϵr = 1+ρₙ*α₂^2*β₆^2*(1+λ)^4/16
    #return F, m, N, ϵr
end

function a_ion(model::MSAIDModel, V, T, z, _data=@f(data))
    σ = model.params.sigma.values
    Z = model.params.charge.values
    nc = length(model)
    isolv1 = findfirst(iszero,Z)
    isolv = isolv1:isolv1
    iions = @iions
    (α₀,α₂,Δ,ξ₂,χ,σₙ,ρₙ,ρ,∑z) = _data

    Γ, B, b₂ = @f(solve_MSAID, _data)
    β₆ = 1-b₂/6 # Checked
    λ  = (1+b₂/3)/β₆ # Checked
    y₁ = 4/(β₆*(1+λ)^2) # Checked

    W₁ = ρ*sum(z[i]*Z[i]^2/(β₆*(σₙ+σ[i]*λ)*(1+σ[i]*Γ)) for i ∈ iions)/∑z # Checked
    W₂ = 1/2*ρₙ*ρ*σₙ^2*B*sum(z[i]*σ[i]^2*Z[i]^2/(2*β₆*(σₙ+σ[i]*λ)*(1+σ[i]*Γ))^2 for i ∈ iions)/∑z # Checked
    Vη = (-W₁/2+√((W₁/2)^2+2B*W₂/β₆^2))/(W₂) # Checked

    ΔΓ = @. Vη*ρₙ*σₙ^2*σ^2*B/(8*β₆*(σₙ+λ*σ)) # Checked
    Dᶠ = @. Z*β₆/(2*(1+σ*Γ-ΔΓ)) # Checked

    D = 1 + Vη^2*ρₙ*ρ*σₙ^2*sum(z[i]*σ[i]^2*Dᶠ[i]^2/(2β₆*(σₙ+λ*σ[i]))^2 for i ∈ iions)/∑z # Checked

    Dac = ρ*sum(z[i]*Dᶠ[i]^2 for i ∈ iions)/∑z # Checked
    Ω   = Vη*ρ*sum(z[i]*σ[i]*Dᶠ[i]^2/(σₙ+λ*σ[i]) for i ∈ iions)/∑z # Checked

    Γₛ = ΔΓ
    Γₛ  .= @. ((1+σ*Γ-ΔΓ)*D-1)/σ # Checked
    a¹  = D*β₆*(σₙ*B/2 + Ω*λ/(D*β₆))/(2*Dac) # Checked

    Jp = zero(Base.promote_eltype(model,V,T,z))

    for i ∈ iions
        σi,Zi = σ[i],Z[i]
        for j ∈ iions
            σij = (σi+σ[j])/2
            Q00 = 2π/Δ*(σij+π*σi*σ[j]*ξ₂/(4Δ))
                - 1/2*Dᶠ[i]*Dᶠ[j]*(ρₙ*σₙ^2*Vη^2/(D*β₆^2*(σₙ+λ*σi)*(σₙ+λ*σ[j]))
                + 4*Γₛ[i]*Γₛ[j]/(D*Dac))
            Jp += z[i]*z[j]*σij*Q00^2/(3π)^2
        end
    end

    for i ∈ iions
        σi,Zi = σ[i],Z[i]
        j = isolv1
        σij = (σi+σ[j])/2
        Q00 = 2π/Δ*(σij+π*σi*σ[j]*ξ₂/(4Δ))
        Jp += 2*z[i]*z[j]*σij*Q00^2/(3π)^2

        Q01 = Dᶠ[i]/(D*β₆)*(λ*Vη/(σₙ+λ*σi)+2*Γₛ[i]*a¹)
        Jp += 2/9*z[i]*z[j]*σij*Q01^2/(3π)^2
    end

    for i ∈ isolv
        σi = σ[i]
        for j ∈ isolv
            σij = (σi+σ[j])/2
            Q11 = 2λ/(D*ρₙ*σₙ^2)*(λ+ρₙ*σₙ^2*Ω*a¹/(2*β₆))
                + σₙ*B*a¹/(2*β₆) - 2/(ρₙ*σₙ^2)

            qp = -b₂*(λ+3)/(1+λ)^2/(ρₙ*σₙ^2)

            h11 = (Q11+2qp)/(2*√(3)*π)
            h12 = √(10)*(Q11-qp)/(2*√(3)*π)
            Jp += z[i]*z[j]*σij*(h11^2 + 1/5*h12^2)
        end
    end

    Jp *= ρ/(3π)/∑z/∑z
    ∑xZN = zero(Jp)
    for i ∈ iions
        σi,Zi,zi = σ[i],Z[i],z[i]
        Ni = 2*Dᶠ[i]/(β₆*σi)*(1+Vη*ρₙ*σₙ^3*B*σi/(24*(σₙ+λ*σi))) - Zi/σi
        ∑xZN += zi*Zi*Ni
    end
    ∑xZN /= ∑z
    #N = @. (2*Dᶠ/(β₆*σ)*(1+Vη*ρₙ*σₙ^3*B*σ/(24*(σₙ+λ*σ))) - Z/σ)
    return (α₀^2*∑xZN - ρₙ/ρ*α₀*B*α₂)/(6π) - Jp
end

function solve_MSAID(model::MSAIDModel,V,T,z,_data = @f(data))
    _1 = oneunit(Base.promote_eltype(model,V,T,z))
    Z = model.params.charge.values
    (α₀,α₂,Δ,ξ₂,χ,σₙ,ρₙ,ρ,∑z) = _data
    Γ0 = sqrt(4π*e_c^2/(4π*ϵ_0*80*k_B*T)*ρ*sum(z[i]*Z[i]^2 for i ∈ @iions)/∑z)*1e-9
    x0 = SVector((Γ0,_1*1.0,_1*2.02)) #TOOD: any better initial point?
    f(x) = obj_MSAID(model,z,x[1]*1e9,x[2]*1e17,x[3],_data)
    _x = Solvers.nlsolve2(f, x0, Solvers.Newton2Var())
    #_x = Solvers.x_sol(sol)
    Γ, B, b₂ = _x[1]*1e9, _x[2]*1e17, _x[3]
    return Γ, B, b₂
end

function dielectric_constant(model::MSAIDModel, V, T, z, _data = @f(data))
    (α₀,α₂,Δ,ξ₂,χ,σₙ,ρₙ,ρ,x) = _data
    Γ, B, b₂ = @f(solve_MSAID, _data)
    β₆ = 1-b₂/6 # Checked
    λ  = (1+b₂/3)/β₆ # Checked
    ϵr = 1 + ρₙ*α₂^2 * β₆^2*(1 + λ)^4 / 16
end