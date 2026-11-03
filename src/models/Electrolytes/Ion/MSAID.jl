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
        charge = nothing,
        RSPmodel = nothing,
        userlocations = String[],
        RSPmodel_userlocations = nothing,
        verbose = false)

## Input parameters

- `sigma`: Single Parameter (`Float64`) — Hard-sphere diameter `[m]`
- `dipole`: Single Parameter (`Float64`) — Permanent dipole moment of the solvent `[C·m]`
- `charge`: Single Parameter (`Float64`) — Charge number `[-]`

## Description 

Mean Spherical Approximation (MSA) ion–dipole model for electrolyte solutions.

The MSAID term adds the electrostatic contribution to the residual Helmholtz energy of a mixture of charged hard spheres (ions) and hard dipoles (solvent), all interacting through the MSA closure. 
The model is built from two species types:

- one neutral solvent with a permanent point dipole `dipole` and hard-sphere diameter `sigma`,
- any number of ions with hard-sphere diameter `sigma` and charge number `charge` (in units of the elementary charge).

Three internal parameters carry the whole electrostatic response:

- `Γ`  — ionic screening parameter (MSA analogue of the Debye `κ`),
- `B`  — ion–dipole coupling parameter,
- `b₂` — dipole–dipole coupling parameter.

They are obtained by solving a closed system of three algebraic equations (Blum 1974 for equal diameters; Blum & Wei 1987 for the general-size case).
All thermodynamic properties are then expressed as closed-form functions of `(Γ, B, b₂)`.

## Model construction

- Exactly one neutral solvent is required; the model errors otherwise.
- `dipole` is only read for the neutral component.
- If a `b` (covolume) parameter is present instead of `sigma`, the hard-sphere diameter is recovered from `b = (3/2) N_A σ³` before the MSA equations are built.

!!! note "unequal ion diameters"
    The Blum–Wei (1987) general-size solution requires solving a coupled inner system before the outer closures in `(Γ, B, b₂)` can be evaluated. 
    This implementation replaces that inner system with an approximation that behaves correctly on the dilute limit.
    The approximation implemented is exact for equal size ion diameters; accuracy is expected to be good for small size spreads, degrading as the ionic diameter spread widens.

## References

1. Blum, L. (1974). Solution of a model for the solvent-electrolyte interactions in the mean spherical approximation. The Journal of Chemical Physics, 61(5), 2129–2133. [doi:10.1063/1.1682224](https://doi.org/10.1063/1.1682224)
1. Blum, L., & Wei, D. Q. (1987). Analytical solution of the mean spherical approximation for an arbitrary mixture of ions in a dipolar solvent. The Journal of Chemical Physics, 87(1), 555–565. [doi:10.1063/1.453604](https://doi.org/10.1063/1.453604)
"""
function MSAID(solvents,ions; charge = nothing, userlocations, verbose=false, RSPmodel = nothing, RSPmodel_userlocations = nothing)
    solvents = format_components(solvents)
    ions = format_components(ions)
    components = vcat(solvents, ions)
    
    userlocations = normalize_userlocations(userlocations)

    params = getparams(components, String[]; userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    if any(keys(params).=="b")
        params["b"].values .*= 3/2/N_A/π*1e-3
        params["b"].values .^= 1/3
        sigma = SingleParam("sigma",components,params["b"].values)
    else
        params["sigma"].values .*= 1E-10
        sigma = params["sigma"]
    end

    dipole = params["dipole"]
    dipole.values .*= 1/(𝕔)*1e-21

    Z = init_charge(components,charge;verbose)

    packagedparams = MSAIDParam(sigma,dipole,SingleParam("charge",components,Z))

    references = String["10.1063/1.1682224"]
    nsolv = count(iszero,Z)
    if !isone(nsolv)
        throw(error("MSAID only supports one neutral solvent."))
    end
    model = MSAID(components, packagedparams, references)
    return model
end

IonDependency(model::MSAIDModel) = IndependentIonModel()

function data(model::MSAIDModel, V, T, z , iondata = (model.params.charge.values, model.params.sigma.values, 1.0))
    β = 1/(k_B*T)
    σ = model.params.sigma.values
    #Z = model.params.charge.values
    Z, _, _ = iondata
    isolv = findfirst(iszero,Z)::Int
    μ = model.params.dipole.values[isolv]
    ∑z = sum(z)
    ρ = N_A*∑z/V
    σₙ = σ[isolv]
    xₙ = z[isolv]/∑z
    ρₙ = ρ*xₙ
    α₀ = e_c*sqrt(β/ϵ_0) # Checked
    α₂ = μ*sqrt(β/3/ϵ_0) # Checked
    _data = (α₀,α₂,σₙ,ρₙ,ρ,σ,Z,z,∑z)
    return _data
end

obj_MSAID(x,_data) = obj_MSAID(x[1],x[2],x[3],_data)

function obj_MSAID(Γ,B,b₂,_data)
    
    (α₀,α₂,σₙ,ρₙ,ρ,σ,Z,z,∑z) = _data
    iions = iions_iterator(Z)

    β₆ = 1-b₂/6 # Checked
    λ  = (1+b₂/3)/β₆ # Checked
    y₁ = 4/(β₆*(1+λ)^2) # Checked

    W₁ = ρ*sum(z[i]*Z[i]^2/(β₆*(σₙ + σ[i]*λ)*(1+σ[i]*Γ)) for i ∈ iions)/∑z # Checked
    W₂ = 1/2*ρₙ*ρ*σₙ^2*B*sum(z[i]*σ[i]^2*Z[i]^2/(2*β₆*(σₙ+σ[i]*λ)*(1+σ[i]*Γ))^2 for i ∈ iions)/∑z # Checked
    Vη = (-W₁/2+sqrt((W₁/2)^2+2B*W₂/β₆^2))/(W₂) # Checked


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


function a_res(model::MSAIDModel, V, T, z, iondata = (model.params.charge.values, model.params.sigma.values) , _data = @f(data,iondata))
    (α₀,α₂,σₙ,ρₙ,ρ,σ,Z,z,∑z) = _data

    isolv1 = findfirst(iszero,Z)
    nc = length(model)
    isolv = isolv1:isolv1    
    iions = iions_iterator(Z)

    Γ, B, b₂ = @f(solve_MSAID, _data)

    Δ = 1 - π*ρ/6*sum(z[i]*σ[i]^3 for i ∈ 1:nc)/∑z
    ξ₂ = ρ*sum(z[i]*σ[i]^2 for i ∈ 1:nc)/∑z
    β₆ = 1-b₂/6 # Checked
    λ  = (1+b₂/3)/β₆ # Checked
    #y₁ = 4/(β₆*(1+λ)^2) # Checked

    W₁ = ρ*sum(z[i]*Z[i]^2/(β₆*(σₙ+σ[i]*λ)*(1+σ[i]*Γ)) for i ∈ iions)/∑z # Checked
    W₂ = 1/2*ρₙ*ρ*σₙ^2*B*sum(z[i]*σ[i]^2*Z[i]^2/(2*β₆*(σₙ+σ[i]*λ)*(1+σ[i]*Γ))^2 for i ∈ iions)/∑z # Checked
    Vη = (-W₁/2+sqrt((W₁/2)^2+2B*W₂/β₆^2))/(W₂) # Checked

    ΔΓ = @. Vη*ρₙ*σₙ^2*σ^2*B/(8*β₆*(σₙ+λ*σ)) # Checked
    Dᶠ = @. Z*β₆/(2*(1+σ*Γ-ΔΓ)) # Checked

    D = 1 + Vη^2*ρₙ*ρ*σₙ^2*sum(z[i]*σ[i]^2*Dᶠ[i]^2/(2β₆*(σₙ+λ*σ[i]))^2 for i ∈ iions)/∑z # Checked

    Dac = ρ*sum(z[i]*Dᶠ[i]^2 for i ∈ iions)/∑z # Checked
    Ω   = Vη*ρ*sum(z[i]*σ[i]*Dᶠ[i]^2/(σₙ+λ*σ[i]) for i ∈ iions)/∑z # Checked

    Γₛ = ΔΓ # in-place reuse of ΔΓ buffer; ΔΓ not used after this
    Γₛ  .= @. ((1+σ*Γ-ΔΓ)*D-1)/σ # Checked
    a¹  = D*β₆*(σₙ*B/2 + Ω*λ/(D*β₆))/(2*Dac) # Checked

    Jp = zero(Base.promote_eltype(model,V,T,z))

    for i ∈ iions
        σi,zi,Dᶠi,Γₛi = σ[i],z[i],Dᶠ[i],Γₛ[i]
        for j ∈ iions
            σij = (σi+σ[j])/2
            Q00 = 2π/Δ*(σij+π*σi*σ[j]*ξ₂/(4Δ)) - 1/2*Dᶠi*Dᶠ[j]*(ρₙ*σₙ^2*Vη^2/(D*β₆^2*(σₙ+λ*σi)*(σₙ+λ*σ[j])) + 4*Γₛi*Γₛ[j]/(D*Dac))
            Jp += zi*z[j]*σij*Q00^2/(3π)^2
        end
    end

    for i ∈ iions
        σi = σ[i]
        j = isolv1
        σij = (σi+σ[j])/2
        Q00 = 2π/Δ*(σij+π*σi*σ[j]*ξ₂/(4Δ))
        Jp += 2*z[i]*z[j]*σij*Q00^2/(3π)^2
        Q01 = Dᶠ[i]/(D*β₆)*(λ*Vη/(σₙ+λ*σi)+2*Γₛ[i]*a¹)
        Jp += 2/9*z[i]*z[j]*σij*Q01^2/(3π)^2
    end

    for i ∈ isolv
        σi,zi = σ[i],z[i]
        for j ∈ isolv
            σij = (σi+σ[j])/2
            #this was bugged before: extra line was independent expression
            Q11 = 2λ/(D*ρₙ*σₙ^2)*(λ+ρₙ*σₙ^2*Ω*a¹/(2*β₆)) + σₙ*B*a¹/(2*β₆) - 2/(ρₙ*σₙ^2)

            qp = -b₂*(λ+3)/(1+λ)^2/(ρₙ*σₙ^2)

            h11 = (Q11+2qp)/(2*√(3)*π)
            h12 = √(10)*(Q11-qp)/(2*√(3)*π)
            Jp += zi*z[j]*σij*(h11^2 + 1/5*h12^2)
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

solve_MSAID(model::MSAIDModel,V,T,z,_data) = solve_MSAID(_data)
solve_MSAID(model::MSAIDModel,V,T,z) = solve_MSAID(data(model,V,T,z))

function solve_MSAID(_data)
    if any(has_dual,_data)
        λdata = primalval(_data)
        λresult = solve_MSAID_impl(λdata)
        result  = solve_MSAID_ad(λresult,_data,λdata)
    else
        result = solve_MSAID_impl(_data)
    end
    return result
end

function solve_MSAID_ad(result::RES,∂data::TUP1,λdata::TUP2) where {RES,TUP1,TUP2}
    if any(has_dual,∂data) # do check here to avoid recomputation if no AD
        ∂result = __gradients_for_root_finders(result,∂data,λdata,obj_MSAID)
        return ∂result
    end
    return result
end

function solve_MSAID_impl(_data)
    (α₀,α₂,σₙ,ρₙ,ρ,σ,Z,z,∑z) = _data
    iions = iions_iterator(Z)

    # ------------------------------------------------------------------
    # (1) b₂ from pure-dipole MSA (Blum 1974, Eq. 3.7)
    # ------------------------------------------------------------------
    d₂  = sqrt(ρₙ) * α₂
    b₂0 = 12 * d₂^2 / (24 + 4 * d₂^2)          # ∈ (0, 3)

    # ------------------------------------------------------------------
    # (2) ε_w from Wertheim MSA at b₂0
    # ------------------------------------------------------------------
    β₆0  = 1 - b₂0/6
    β₃0  = 1 + b₂0/3
    β₁₂0 = 1 + b₂0/12
    λ0   = β₃0/β₆0
    εw0  = 1 + ρₙ * α₂^2 * β₆0^2 * (1 + λ0)^4 / 16

    # ------------------------------------------------------------------
    # (3) Γ* from Debye — no T, no k_B, no ϵ_0, no e_c
    #     Γ0_si² = α₀² · ρ · Σ(zᵢ Zᵢ²)/Σz / εw0
    # ------------------------------------------------------------------
    I_ion   = sum(z[i]*Z[i]^2 for i ∈ iions) / ∑z
    Γ0_si   = α₀ * sqrt( ρ * I_ion / εw0 )
    Γ0_star = Γ0_si * σₙ

    # ------------------------------------------------------------------
    # (4) B* from dilute Blum-Wei expansion
    #     All factors already in _data; no explicit T either.
    # ------------------------------------------------------------------
    σ̄ = sum(z[i]*σ[i] for i ∈ iions) / sum(z[i] for i ∈ iions)
    B0_si = Γ0_si^2 * sqrt(εw0) * (α₂/α₀) * β₁₂0^2 * β₆0 /
            ( (β₆0 + σ̄ * β₃0 / σₙ) * σₙ )
    B0_star = B0_si * σₙ^2

    # ------------------------------------------------------------------
    # (5) dimensionless initial point
    # ------------------------------------------------------------------
    x0 = SVector(promote(Γ0_star, B0_star, b₂0))
    TT = eltype(x0)

    # ------------------------------------------------------------------
    # (6) residual in dimensionless coordinates
    # ------------------------------------------------------------------
    f(x) = obj_MSAID(x[1] / σₙ,          # Γ SI
                     x[2] / σₙ^2,        # B SI
                     x[3],               # b₂
                     _data)

    # ------------------------------------------------------------------
    # (7) bounds
    # ------------------------------------------------------------------
    lb = (zero(TT), zero(TT), zero(TT))
    ub = (TT(1.0e4), TT(1.0e3), TT(3.0))

    _x = Solvers.nlsolve2(f, x0, Solvers.Newton2Var();
                          bounds = (lb, ub))

    Γ  = _x[1] / σₙ       # SI, m⁻¹
    B  = _x[2] / σₙ^2     # SI, m⁻²
    b₂ = _x[3]
    return SVector((Γ, B, b₂))
end
function dielectric_constant(model::MSAIDModel, V, T, z, _data = @f(data))
    (α₀,α₂,σₙ,ρₙ,ρ,σ,Z,z,∑z) = _data
    Γ, B, b₂ = @f(solve_MSAID, _data)
    β₆ = 1-b₂/6 # Checked
    λ  = (1+b₂/3)/β₆ # Checked
    ϵr = 1 + ρₙ*α₂^2 * β₆^2*(1 + λ)^4 / 16
    return ϵr
end