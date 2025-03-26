struct BornParam <: EoSParam
    sigma_born::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type BornModel <: EoSModel end

struct Born{ϵ} <: BornModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::BornParam
    RSPmodel::ϵ
    references::Array{String,1}
end

export Born

"""
    Born(solvents::Array{String,1},
        salts::Array{String,1};
        RSPmodel = ConstRSP,
        userlocations = String[],
        RSPmodel_userlocations = String[],
        verbose = false)

## Input parameters
- `sigma_born`: Single Parameter (`Float64`) - Born Diameter `[m]`
- `charge`: Single Parameter (`Float64`) - Charge `[-]`

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a Born model. The Born term gives the excess Helmholtz energy to account for the electrostatic interactions between ions in solution.

## References
1. Born, M. (1920). Z. Phys. 1, 45.
"""
function Born(solvents,salts; RSPmodel=ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)
    _solvents = group_components(solvents)
    ions = ion_groups.flattenedgroups
    components = deepcopy(ions)
    prepend!(components,_solvents)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    params = getparams(components, ["Electrolytes/Born/Born.csv","Electrolytes/properties/charges.csv"]; userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    params["sigma_born"].values .*= 1E-10
    sigma_born = params["sigma_born"]
    charge = params["charge"]

    packagedparams = BornParam(sigma_born,charge)

    references = String[]

    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    model = Born(components, _solvents, ions, isolvents, iions, packagedparams, init_RSPmodel,references)
    return model
end

a_res(::Nothing,V,T,z,_data=nothing) = 0.0

function data(model::BornModel, V, T, z)
    return dielectric_constant(model, V, T, z)
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