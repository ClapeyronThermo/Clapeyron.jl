abstract type WAvgILModel <: RSPModel end

struct WAvgILParam <: EoSParam
    e_r::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct WAvgIL <: WAvgILModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::WAvgILParam
    references::Array{String,1}
end

@registermodel WAvgIL
export WAvgIL
function WAvgIL(solvents,salts; userlocations::Vector{String}=String[],assoc_userlocations::Vector{String}=String[], verbose::Bool=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,[salts[1][1]])
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    
    params = getparams(components, ["Electrolytes/RSP/WAvg_like.csv"]; userlocations=userlocations, verbose=verbose)
    e_r = params["e_r"]
    Mw  = params["Mw"]
    packagedparams = WAvgILParam(e_r,Mw)

    references = String[]
    
    model = WAvgIL(components, solvents, ions, isolvents, iions, packagedparams,references)
    return model
end

function dielectric_constant(model::WAvgILModel, V, T, z,_data = nothing)
    #z_salt = ∑(z[model.iions])./length(model.iions)
    ϵᵣ = model.params.e_r.values
    Mw = model.params.Mw.values

    ninv = 1/length(model.iions)
    ∑z = sum(z)
    ∑zinv = 1/∑z
    res = zero(eltype(z))
    ∑wmw = res
    w = res
    for i ∈ model.isolvents
        w = z[i]*Mw[i]*∑zinv
        ∑wmw += w
        res += w*ϵᵣ[i]
    end
    
    for i ∈ model.iions
        w = z[i]*Mw[i]*∑zinv*ninv
        ∑wmw += w
        res += w*ϵᵣ[i]
    end

    return res/∑wmw
    #z = append!(z[model.isolvents],z_salt)
    #x = z./∑(z)

    #w = x.*Mw./∑(x.*Mw);

    #return ∑(w.*ϵᵣ)
end

is_splittable(::WAvgIL) = false