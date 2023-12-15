struct DHParam <: EoSParam
    sigma::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type DHModel <: IonModel end

struct DH{ϵ} <: DHModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::DHParam
    RSPmodel::ϵ
    references::Array{String,1}
end

@registermodel DH

"""
    DH(solvents,salts; 
    RSPmodel=ConstW, 
    SAFTlocations=String[], 
    userlocations=String[], 
    ideal_userlocations=String[], 
    verbose=false)

Debye-Hückel (DH) model for electrostatic interaction.
"""
DH

export DH
function DH(solvents,ions; RSPmodel=ConstW, SAFTlocations=String[], userlocations=String[], verbose=false)

    components = deepcopy(ions)
    prepend!(components,solvents)
    icomponents = 1:length(components)

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

    init_RSPmodel = RSPmodel(solvents,ions)

    model = DH(components, icomponents, packagedparams, init_RSPmodel,references)
    return model
end

function data(model::DHModel, V, T, z)
    return dielectric_constant(model.RSPmodel, V, T, z)
end

function a_res(model::DHModel, V, T, z,_data=@f(data))  
    ϵ_r = _data
    σ = model.params.sigma.values
    Z = model.params.charge.values
    iions = model.icomponents[Z.!=0]

    if length(iions) == 0
        return zero(V+T+first(z))
    end
   
    ∑z = sum(z)
    #x = z ./ sum(z)
    ρ = N_A*sum(z)/V

    s = e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)
    κ = (4π*s*ρ*sum(z[i]*Z[i]^2 for i ∈ iions)/∑z)^(1/2)
    y = σ*κ
    χ = @. 3/y^3*(3/2+log1p(y)-2*(1+y)+1/2*(1+y)^2)
    return -1/3*s*κ*sum(z[i]*Z[i]^2*χ[i] for i ∈ iions)/∑z
end