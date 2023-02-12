struct BornParam <: EoSParam
    sigma_born::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type BornModel <: EoSModel end

struct Born{ϵ} <: BornModel
    components::Array{String,1}
    solvents::Union{Array{String,1},Array{Any,1}}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::BornParam
    RSPmodel::ϵ
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel Born
export Born
function Born(solvents,salts; RSPmodel=ConstW, SAFTlocations=String[], userlocations=String[], ideal_userlocations=String[], verbose=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(ions)
    prepend!(components,solvents)
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)
    
    params,sites = getparams(components, append!(["Electrolytes/Born/Born.csv","Electrolytes/properties/charges.csv","properties/molarmass.csv"],SAFTlocations); userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    icomponents = 1:length(components)
    params["sigma_born"].values .*= 1E-10
    sigma_born = params["sigma_born"]
    charge = params["charge"]

    packagedparams = BornParam(sigma_born,charge)

    references = String[]

    if RSPmodel !== nothing
        init_RSPmodel = RSPmodel(solvents,salts)
    else
        init_RSPmodel = nothing
    end

    model = Born(components, solvents, ions, icomponents, isolvents, iions, packagedparams, init_RSPmodel, 1e-12,references)
    return model
end

a_res(::Nothing,V,T,z,_data=nothing) = 0.0

function data(model::BornModel, V, T, z)
    return dielectric_constant(model.RSPmodel, V, T, z)
end

function a_res(model::BornModel, V, T, z,_data=@f(data))
    σ_born = model.params.sigma_born.values
    if length(model.iions) == 0
        return zero(T+first(z))
    end
    Z = model.params.charge.values
    ϵ_r = _data
    return -e_c^2/(4π*ϵ_0*k_B*T*sum(z))*(1-1/ϵ_r)*sum(z[i]*Z[i]^2/σ_born[i] for i ∈ model.iions)
end