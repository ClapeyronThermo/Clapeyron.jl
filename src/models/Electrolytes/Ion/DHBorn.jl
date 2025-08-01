

abstract type DHBornModel <: DHModel end

struct DHBorn{ϵ} <: DHBornModel
    components::Array{String,1}
    params::BornParam
    RSPmodel::ϵ
    references::Array{String,1}
end

export DHBorn
"""
    DHBorn(solvents::Array{String,1},
        ions::Array{String,1};
        RSPmodel = ConstRSP,
        userlocations = String[],
        RSPmodel_userlocations = String[],
        verbose = false)

## Input parameters
- `sigma`: Single Parameter (`Float64`) - Diameter of closest approach `[m]`
- `sigma_born`: Single Parameter (`Float64`) - Born Diameter `[m]`
- `charge`: Single Parameter (`Float64`) - Charge `[-]`

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a Debye-Hückel-Born model. The Debye-Hückel-Born term gives the excess Helmholtz energy to account for the electrostatic interactions between ions in solution.

## References
1. Debye, P., Huckel, E. (1923). Phys. Z. 24, 185.
2. Born, M. (1920). Z. Phys. 1, 45.
"""
function DHBorn(solvents,ions; RSPmodel=ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)
    components = deepcopy(ions)
    prepend!(components,solvents)
    params = getparams(components, append!(["Electrolytes/Born/born_like.csv"]); userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    sigma_born = params["sigma_born"]
    sigma_born.values .*= 1E-10

    packagedparams = BornParam(sigma_born)

    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    references = String[]

    model = DHBorn(components, packagedparams, init_RSPmodel, references)
    return model
end

function a_res(model::DHBornModel, V, T, z, iondata)
    return a_dh(model, V, T, z, iondata) + a_born(model, V, T, z, iondata)
end
