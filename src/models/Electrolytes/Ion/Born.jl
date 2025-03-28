struct BornParam <: EoSParam
    sigma_born::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type BornModel <: EoSModel end

struct Born{ϵ} <: BornModel
    components::Array{String,1}
    params::BornParam
    RSPmodel::ϵ
    references::Array{String,1}
end

export Born

"""
    Born(solvents::Array{String,1},
        ions::Array{String,1};
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
function Born(solvents,ions; RSPmodel = ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)
    
    components = deepcopy(ions)
    prepend!(components,solvents)    
    params = getparams(components, ["Electrolytes/Born/Born.csv","Electrolytes/properties/charges.csv"]; userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    params["sigma_born"].values .*= 1E-10
    sigma_born = params["sigma_born"]
    charge = params["charge"]

    packagedparams = BornParam(sigma_born,charge)

    references = String[]

    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    model = Born(components, packagedparams, init_RSPmodel, references)
    return model
end


function data(model::BornModel, V, T, z)
    return dielectric_constant(model, V, T, z)
end

function a_res(model::BornModel,V,T,z, _data = @f(data))
    return a_born(model,V,T,z,_data)
end

function a_born(model::IonModel, V, T, z, ϵ_r = dielectric_constant(model,V,T,z),σ_born = model.params.sigma_born.values)
    Z = model.params.charge.values
    ϵ_r = _data
    s = -e_c^2/(4π*ϵ_0*k_B*T*sum(z))
    if all(iszero,Z)
        return zero(Base.promote_eltype(s,σ_born))
    end
    return s*(1-1/ϵ_r)*sum(z[i]*Z[i]*Z[i]/σ_born[i] for i ∈ @iions)
end