struct DHParam <: EoSParam
    sigma::SingleParam{Float64}
    charge::SingleParam{Float64}
end

abstract type DHModel <: IonModel end

struct DH{ϵ} <: DHModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
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
function DH(solvents,salts; RSPmodel=ConstW, SAFTlocations=String[], userlocations=String[], ideal_userlocations=String[], verbose=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(ions)
    _solvents = group_components(solvents)
    prepend!(components,_solvents)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    if occursin("CPA",SAFTlocations[1])
        params,sites = getparams(components, append!(["Electrolytes/properties/charges.csv","properties/molarmass.csv"],SAFTlocations); userlocations=userlocations,ignore_missing_singleparams=["b","charge"], verbose=verbose)
        params["b"].values .*= 3/2/N_A/π*1e-3
        params["b"].values .^= 1/3
        sigma = params["b"]
    else
        params,sites = getparams(components, append!(["Electrolytes/properties/charges.csv","properties/molarmass.csv"],SAFTlocations); userlocations=userlocations,ignore_missing_singleparams=["sigma","charge"], verbose=verbose)
        params["sigma"].values .*= 1E-10
        sigma = params["sigma"]
    end
    
    charge = params["charge"]

    packagedparams = DHParam(sigma,charge)

    references = String[]

    if RSPmodel !== nothing
        init_RSPmodel = RSPmodel(solvents,salts)
    else
        init_RSPmodel = nothing
    end

    model = DH(components, _solvents, ions, isolvents, iions, packagedparams, init_RSPmodel,references)
    return model
end

function data(model::DHModel, V, T, z)
    return dielectric_constant(model.RSPmodel, V, T, z)
end

function a_res(model::DHModel, V, T, z,_data=@f(data))  
    if length(model.iions) == 0
        return zero(V+T+first(z))
    end
    ϵ_r = _data
    σ = model.params.sigma.values
    Z = model.params.charge.values
    ∑z = sum(z)
    #x = z ./ sum(z)
    ρ = N_A*sum(z)/V

    s = e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)
    κ = (4π*s*ρ*sum(z[i]*Z[i]^2 for i ∈ model.iions)/∑z)^(1/2)
    y = σ*κ
    χ = @. 3/y^3*(3/2+log1p(y)-2*(1+y)+1/2*(1+y)^2)
    return -1/3*s*κ*sum(z[i]*Z[i]^2*χ[i] for i ∈ model.iions)/∑z
end