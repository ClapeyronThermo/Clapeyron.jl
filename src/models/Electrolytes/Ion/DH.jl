struct DHParam <: EoSParam
    sigma::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type DHModel <: IonModel end

struct DH{ϵ} <: DHModel
    components::Array{String,1}
    params::DHParam
    RSPmodel::ϵ
    references::Array{String,1}
end

"""
    DH(solvents::Array{String,1},
        ions::Array{String,1};
        RSPmodel = ConstRSP,
        userlocations = String[],
        RSPmodel_userlocations = String[],
        verbose = false)

## Input parameters
- `sigma`: Single Parameter (`Float64`) - Diameter of closest approach `[m]`
- `charge`: Single Parameter (`Float64`) - Charge `[-]`

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a Debye-Hückel model. The Debye-Hückel term gives the excess Helmholtz energy to account for the electrostatic interactions between ions in solution.

## References
1. Debye, P., Huckel, E. (1923). Phys. Z. 24, 185.
"""
DH

export DH
function DH(solvents,ions; RSPmodel=ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)

    components = deepcopy(ions)
    prepend!(components,solvents)

    params = getparams(components, ["Electrolytes/properties/charges.csv","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)

    if any(keys(params).=="b")
        params["b"].values .*= 3/2/N_A/π*1e-3
        params["b"].values .^= 1/3
        sigma = SingleParam("sigma",components,params["b"].values)
    else
        params["sigma"].values .*= 1E-10
        sigma = params["sigma"]
    end

    charge = params["charge"]

    packagedparams = DHParam(sigma,charge)

    references = String[]

    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    model = DH(components, packagedparams, init_RSPmodel,references)
    return model
end

function data(model::DHModel, V, T, z)
    return dielectric_constant(model, V, T, z), model.params.sigma.values
end

function a_res(model::DHModel, V, T, z, _data=@f(data))
    return a_ion(model, V, T, z, _data)
end

function a_ion(model::DHModel, V, T, z,_data=@f(data))
    ϵ_r, σ = _data
    Z = model.params.charge.values

    if all(iszero,Z)
        return zero(V+T+first(z))
    end
    nc = length(model)
    ∑z = sum(z)
    ρ = N_A*sum(z)/V
    s = e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)
    κ = sqrt(4π*s*ρ*sum(z[i]*Z[i]*Z[i] for i ∈ 1:nc)/∑z)
    res = zero(Base.promote_eltype(model,V,T,z))
    for i ∈ @iions
        yi = σ[i]*κ
        yip1 = yi + 1
        χi = 3/(yi*yi*yi)*(3/2+log1p(yi)-2*yip1+1/2*yip1*yip1)
        res +=z[i]*Zi*Zi*χi
    end

    return -1/3*s*κ*res/∑z
    #y = σ*κ
    #χ = @. 3/y^3*(3/2+log1p(y)-2*(1+y)+1/2*(1+y)^2)
    # return -1/3*s*κ*sum(z[i]*Z[i]^2*χ[i] for i ∈ iions)/∑z
end