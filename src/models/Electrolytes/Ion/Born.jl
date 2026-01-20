struct BornParam <: EoSParam
    sigma_born::SingleParam{Float64}
end

abstract type BornModel <: IonModel end

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

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a Born model. The Born term gives the excess Helmholtz energy to account for the electrostatic interactions between ions in solution.

## References
1. Born, M. (1920). Z. Phys. 1, 45.
"""
function Born(solvents,ions; RSPmodel = ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)

    solvents = format_components(solvents)
    ions = format_components(ions)
    components = vcat(solvents, ions)

    userlocations = normalize_userlocations(userlocations)
    RSPmodel_userlocations = normalize_userlocations(RSPmodel_userlocations)

    params = getparams(components, ["Electrolytes/Born/born_like.csv"]; userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    params["sigma_born"].values .*= 1E-10
    sigma_born = params["sigma_born"]

    packagedparams = BornParam(sigma_born)

    references = String[]

    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    model = Born(components, packagedparams, init_RSPmodel, references)
    return model
end


function data(model::BornModel, V, T, z)
    return dielectric_constant(model, V, T, z)
end

function a_res(ionmodel::BornModel, V, T, z, iondata)
    return a_born(ionmodel, V, T, z, iondata)
end

function a_born(model::EoSModel, V, T, z, iondata)
    σ_born = model.params.sigma_born.values
    return a_born(V, T, z, iondata, σ_born)
end

function a_born(V::Number, T, z, iondata, σ_born)
    Z, σ, ϵ_r = iondata
    s = -e_c^2/(4π*ϵ_0*k_B*T*sum(z))
    if all(iszero,Z)
        return zero(Base.promote_eltype(V, T, z, Z, σ, ϵ_r))
    end
    return s*(1-1/ϵ_r)*sum(z[i]*Z[i]*Z[i]/σ_born[i] for i ∈ @iions)
end