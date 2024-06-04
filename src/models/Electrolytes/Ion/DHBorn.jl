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

"""
    DHBorn(solvents,ions;
    RSPmodel=ConstRSP,
    userlocations=String[],
    SAFT_userlocations=String[],
    RSP_userlocations = String[]
    verbose=false)

Debye-Hückel (DH) model for electrostatic interaction with Born term.

Requires `sigma`, that should be provided by an EoS.
"""


export DHBorn
function DHBorn(solvents,ions; RSPmodel=ConstRSP, userlocations=String[], RSP_userlocations=String[], verbose=false)
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

    init_RSPmodel = RSPmodel(solvents,ions)

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
    iions = model.icomponents[Z.!=0]

    if length(iions) == 0
        return zero(T+first(z))
    end
    
    ϵ_r, σ, σ_born = _data

    return -e_c^2/(4π*ϵ_0*k_B*T*sum(z))*(1-1/ϵ_r)*sum(z[i]*Z[i]^2/σ_born[i] for i ∈ iions)
end