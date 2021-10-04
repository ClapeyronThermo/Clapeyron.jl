struct DHParam <: EoSParam
    sigma::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type IonModel <: EoSModel end
abstract type DHModel <: IonModel end

struct DH{ϵ<:RSPModel} <: DHModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::DHParam
    RSPmodel::ϵ
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel DH
export DH
function DH(solvents,salts; RSPmodel=ConstW, SAFTlocations=String[], userlocations=String[], ideal_userlocations=String[], verbose=false)
    ion_groups = GroupParam(salts, ["Electrolytes/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    if occursin("CPA",SAFTlocations[1])
        params,sites = getparams(components, append!(["Electrolytes/charges.csv","properties/molarmass.csv"],SAFTlocations); userlocations=userlocations,ignore_missing_singleparams=["b","charge"], verbose=verbose)
        params["b"].values .*= 3/2/N_A/π*1e-3
        params["b"].values .^= 1/3
        sigma = params["b"]
    else
        params,sites = getparams(components, append!(["Electrolytes/charges.csv","properties/molarmass.csv"],SAFTlocations); userlocations=userlocations,ignore_missing_singleparams=["sigma","charge"], verbose=verbose)
        params["sigma"].values .*= 1E-10
        sigma = params["sigma"]
    end
    
    charge = params["charge"]

    packagedparams = DHParam(sigma,charge)

    references = [""]

    init_RSPmodel = RSPmodel(solvents,salts)

    model = DH(components, solvents, ions, icomponents, isolvents, iions, packagedparams, init_RSPmodel, 1e-12,references)
    return model
end

function a_ion(model::DHModel, V, T, z)
    σ = model.params.sigma.values
    Z = model.params.charge.values
    ϵ_r = RSP(model.RSPmodel,V,T,z)

    x = z ./ sum(z)
    ρ = N_A*sum(z)/V

    s = e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)
    κ = (4π*s*ρ*sum(x[i]*Z[i]^2 for i ∈ model.iions))^(1/2)
    y = σ*κ
    χ = @. 3/y^3*(3/2+log(1+y)-2*(1+y)+1/2*(1+y)^2)
    return -1/3*s*κ*sum(x[i]*Z[i]^2*χ[i] for i ∈ model.iions)
end