abstract type MSABornModel <: MSAModel end

struct MSABorn{ϵ} <: MSABornModel
    components::Array{String,1}
    params::BornParam
    RSPmodel::ϵ
    references::Array{String,1}
end

export MSABorn

"""
    MSABorn(solvents::Array{String,1},
        ions::Array{String,1};
        RSPmodel = ConstRSP,
        userlocations = String[],
        RSPmodel_userlocations = String[],
        verbose = false)

## Input parameters
- `sigma`: Single Parameter (`Float64`) - Hard-sphere diameter `[m]`
- `sigma_born`: Single Parameter (`Float64`) - Born Diameter `[m]`
- `charge`: Single Parameter (`Float64`) - Charge `[-]`

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a Mean Spherical Approximation-Born model. The MSA-Born term gives the excess Helmholtz free energy to account for the electrostatic interactions between ions in solution.

## References
1. Blum, L. (1974). Solution of a model for the solvent‐electrolyte interactions in the mean spherical approximation, 61, 2129–2133.
2. Born, M. (1920). Z. Phys. 1, 45.
"""
function MSABorn(solvents,ions; RSPmodel=ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)
    params = getparams(components, ["Electrolytes/properties/charges.csv","Electrolytes/Born/born_like.csv"]; userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    sigma_born = params["sigma_born"]
    sigma_born.values .*= 1E-10

    packagedparams = BornParam(sigma_born)

    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    references = String[]

    model = MSABorn(components, packagedparams, init_RSPmodel, references)
    return model
end

function a_res(model::MSABornModel, V, T, z, iondata)
    return a_MSA(model, V, T, z, iondata) + a_born(model, V, T, z, iondata)
end
