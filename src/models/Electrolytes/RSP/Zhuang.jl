abstract type ZhuangModel <: RSPModel end

struct ZhuangParam <: EoSParam
    mu::SingleParam{Float64}
end

struct Zhuang <: ZhuangModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::ZhuangParam
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel Zhuang
export Zhuang
function Zhuang(solvents,salts; userlocations::Vector{String}=String[],assoc_userlocations::Vector{String}=String[], verbose::Bool=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    params = getparams(solvents, ["Electrolytes/RSP/Zhuang_like.csv"]; userlocations=userlocations, verbose=verbose)
    params["mu"].values .*= 3.33564e-30
    mu = params["mu"]
    packagedparams = ZhuangParam(mu)

    references = [""]
    
    model = Zhuang(components, solvents, ions, icomponents, isolvents, iions, packagedparams, 1e-12,references)
    return model
end

function RSP(model::ZhuangModel, V, T, z)
    μ = model.params.mu.values

    β = 1/(T*k_B)
    ρ = z*N_A/V

    y = sum(β*ρ[i]*μ[i]^2/(3*ϵ_0) for i ∈ model.isolvents)

    B = 3*y*(2y^2+3y+9)/(y^2+6y+9)

    poly = (2,-(1+B),-1)

    return (-poly[2]+sqrt(poly[2]^2-4*poly[1]*poly[3]))/(2*poly[1])
end

is_splittable(::Zhuang) = false