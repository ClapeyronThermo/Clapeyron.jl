struct MSAParam <: EoSParam
    sigma::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type MSAModel <: IonModel end

struct MSA{ϵ<:RSPModel} <: MSAModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::MSAParam
    RSPmodel::ϵ
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel MSA
export MSA
function MSA(solvents,salts; RSPmodel=ConstW, SAFTlocations=String[], userlocations=String[], ideal_userlocations=String[], verbose=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)
    
    params,sites = getparams(components, append!(["Electrolytes/properties/charges.csv","properties/molarmass.csv"],SAFTlocations); userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    icomponents = 1:length(components)
    params["sigma"].values .*= 1E-10
    sigma = params["sigma"]
    charge = params["charge"]

    packagedparams = MSAParam(sigma,charge)

    references = String[]

    init_RSPmodel = RSPmodel(solvents,salts)

    model = MSA(components, solvents, ions, icomponents, isolvents, iions, packagedparams, init_RSPmodel, 1e-12,references)
    return model
end

function a_ion(model::ElectrolyteModel, V, T, z,ion::MSAModel)
    σ = model.params.sigma.values
    Z = model.params.charge.values
    ϵ_r = RSP(model,V,T,z,model.RSPmodel)
    ∑z = sum(z)
    ρ = N_A*sum(z)/V
    Γ = @f(screening_length)
    Δ = 1-π*ρ/6*sum(z[i]*σ[i]^3 for i ∈ model.iions)/∑z
    Ω = 1+π*ρ/(2*Δ)*sum(z[i]*σ[i]^3/(1+Γ*σ[i]) for i ∈ model.iions)/∑z
    Pn = ρ/Ω*sum(z[i]*σ[i]*Z[i]/(1+Γ*σ[i]) for i ∈ model.iions)/∑z

    U_MSA = -e_c^2*V/(4π*ϵ_0*ϵ_r)*(Γ*ρ*sum(z[i]*Z[i]^2/(1+Γ*σ[i]) for i ∈ model.iions)/∑z + π/(2Δ)*Ω*Pn^2)
    return (U_MSA+Γ^3*k_B*T*V/(3π))/(N_A*k_B*T*sum(z))
end

function screening_length(model::MSAModel,V,T,z)
    σ = model.params.sigma.values
    Z = model.params.charge.values
    ϵ_r = RSP(model.RSPmodel,V,T,z)
    ∑z = sum(z)
    #x = z ./ sum(z)
    ρ = N_A*sum(z)/V
    Δ = 1-π*ρ/6*sum(z[i]*σ[i]^3 for i ∈ model.iions)/∑z

    Γold = (4π*e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)*ρ*sum(z[i]*Z[i]^2 for i ∈ model.iions)/∑z)^(1/2)
    _0 = zero(Γold)
    Γnew = _0
    tol  = one(_0)
    iter = 1
    while tol>1e-8 && iter < 100
        Ω = 1+π*ρ/(2*Δ)*sum(z[i]*σ[i]^3/(1+Γold*σ[i]) for i ∈ model.iions)/∑z
        Pn = ρ/Ω*sum(z[i]*σ[i]*Z[i]/(1+Γold*σ[i]) for i ∈ model.iions)/∑z
        #Q = @. (Z-σ^2*Pn*(π/(2Δ)))./(1+Γold*σ)
        ∑Q2x = _0
        for i ∈ model.iions
            Qi = (Z[i]-σ[i]^2*Pn*(π/(2Δ)))/(1+Γold*σ[i])
            ∑Q2x += z[i]*Qi^2
        end
        Γnew = sqrt(π*e_c^2*ρ/(4π*ϵ_0*ϵ_r*k_B*T)*∑Q2x/∑z)
        tol = abs(1-Γnew/Γold)
        Γold = Γnew
        iter += 1
    end
    return Γnew
end