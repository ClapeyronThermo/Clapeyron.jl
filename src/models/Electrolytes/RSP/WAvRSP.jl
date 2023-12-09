abstract type WAvgRSPModel <: RSPModel end

struct WAvgRSPParam <: EoSParam
    charge::SingleParam{Float64}
    dielectric_constant::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct WAvgRSP <: WAvgRSPModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::WAvgRSPParam
    references::Array{String,1}
end

@registermodel WAvgRSP
export WAvgRSP
function WAvgRSP(solvents,ions; userlocations::Vector{String}=String[],assoc_userlocations::Vector{String}=String[], verbose::Bool=false)

    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)
    
    params = getparams(components, ["Electrolytes/RSP/dielectric.csv","Electrolytes/properties/charges.csv","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose, ignore_missing_singleparams=["dielectric"])
    e_r = params["dielectric"]
    Mw  = params["Mw"]
    charge = params["charge"]
    packagedparams = WAvgRSPParam(charge,e_r,Mw)

    references = String[]
    
    model = WAvgRSP(components, icomponents, packagedparams,references)
    return model
end

function dielectric_constant(model::WAvgRSPModel, V, T, z,_data = nothing)
    #z_salt = ∑(z[model.iions])./length(model.iions)
    ϵᵣ = model.params.dielectric_constant.values
    Mw = model.params.Mw.values
    Z = model.params.charge.values

    isolvent = model.icomponents[Z.==0]

    res = zero(eltype(z))
    ∑wmw = res
    w = res
    for i ∈ isolvent
        w = z[i]*Mw[i]
        ∑wmw += w
        res += w*ϵᵣ[i]
    end

    return res/∑wmw
end

is_splittable(::WAvgRSP) = false