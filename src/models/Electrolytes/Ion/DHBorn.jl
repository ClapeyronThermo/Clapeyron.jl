struct DHBornParam <: EoSParam
    sigma::SingleParam{Float64}
    sigma_born::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type DHBornModel <: DHModel end

struct DHBorn{ϵ} <: DHBornModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::DHBornParam
    RSPmodel::ϵ
    references::Array{String,1}
end

@registermodel DHBorn

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
    icomponents = 1:length(components)
    params = getparams(components, append!(["Electrolytes/properties/charges.csv","properties/molarmass.csv","Electrolytes/Born/born_like.csv"]); userlocations=userlocations,ignore_missing_singleparams=["sigma_born","charge"], verbose=verbose)
    if any(keys(params).=="b")
        params["b"].values .*= 3/2/N_A/π*1e-3
        params["b"].values .^= 1/3
        sigma = SingleParam("sigma",components,params["b"].values)
    else
        params["sigma"].values .*= 1E-10
        sigma = params["sigma"]
    end

    charge = params["charge"]
    sigma_born = params["sigma_born"]
    sigma_born.values .*= 1E-10

    packagedparams = DHBornParam(sigma,sigma_born,charge)

    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    references = String[]

    model = DHBorn(components, icomponents, packagedparams, init_RSPmodel, references)
    return model
end

function a_res(model::DHBornModel, V, T, z, _data=@f(data))
    _data_DH = (_data[1], _data[2])
    return a_ion(model, V, T, z, _data_DH)+a_born(model, V, T, z, _data)
end

function data(model::DHBornModel, V, T, z)
    return dielectric_constant(model.RSPmodel, V, T, z), model.params.sigma.values, model.params.sigma_born.values
end

function a_born(model::DHBornModel, V, T, z,_data=@f(data))
    Z = model.params.charge.values
    ϵ_r, σ, σ_born = _data

    if all(iszero,Z)
        return zero(Base.promote_eltype(model,T,z,ϵ_r, σ, σ_born))
    end

    res = zero(Base.promote_eltype(z,Z,σ_born))
    for i in model.icomponents
        if Z[i] != 0
            res += z[i]*Z[i]^2/σ_born[i]
        end
    end

    return -e_c^2/(4π*ϵ_0*k_B*T*sum(z))*(1-1/ϵ_r)*res
end