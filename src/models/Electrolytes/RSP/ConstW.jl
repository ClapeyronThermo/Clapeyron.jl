abstract type ConstWModel <: RSPModel end

struct ConstWParam <: EoSParam
end

struct ConstW <: ConstWModel
    components::Array{String,1}
    solvents::Union{Array{String,1},Array{Any,1}}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    params::ConstWParam
    absolutetolerance::Float64
    references::Array{String,1}
end

@registermodel ConstW
export ConstW
function ConstW(solvents,salts; userlocations::Vector{String}=String[], verbose::Bool=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)
    icomponents = 1:length(components)

    references = String[]
    
    model = ConstW(components, solvents, ions, icomponents, isolvents, iions, ConstWParam(), 1e-12,references)
    return model
end

function dielectric_constant(model::ConstWModel,V,T,z,_data=nothing)
    return 78.4
end

is_splittable(::ConstW) = false