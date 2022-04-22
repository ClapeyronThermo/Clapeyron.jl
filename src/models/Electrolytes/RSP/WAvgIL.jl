abstract type WAvgILModel <: RSPModel end

struct WAvgILParam <: EoSParam
    e_r::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct WAvgIL <: WAvgILModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::WAvgILParam
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel WAvgIL
export WAvgIL
function WAvgIL(solvents,salts; userlocations::Vector{String}=String[],assoc_userlocations::Vector{String}=String[], verbose::Bool=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,[salts[1][1]])
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    
    params = getparams(components, ["Electrolytes/RSP/WAvg_like.csv"]; userlocations=userlocations, verbose=verbose)
    e_r = params["e_r"]
    Mw  = params["Mw"]
    packagedparams = WAvgILParam(e_r,Mw)

    references = [""]
    
    model = WAvgIL(components, solvents, ions, icomponents, isolvents, iions, packagedparams, 1e-12,references)
    return model
end

function RSP(electromodel::ElectrolyteModel, V, T, z,model::WAvgILModel)
    z_salt = ∑(z[model.iions])./length(model.iions)
    z = append!(z[model.isolvents],z_salt)
    x = z./∑(z)
    ϵᵣ = model.params.e_r.values
    Mw = model.params.Mw.values

    w = x.*Mw./∑(x.*Mw);

    return ∑(w.*ϵᵣ)
end

is_splittable(::WAvgIL) = false