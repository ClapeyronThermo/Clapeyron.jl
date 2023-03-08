abstract type ConstWModel <: RSPModel end

struct ConstW <: ConstWModel
    components::Array{String,1}
    solvents::Array{String,1}
    ions::Array{String,1}
    isolvents::UnitRange{Int}
    iions::UnitRange{Int}
    value::Float64
    references::Array{String,1}
end

@registermodel ConstW
export ConstW
function ConstW(solvents,salts; userlocations::Vector{String}=String[], value =  78.4, verbose::Bool=false)
    ion_groups = GroupParam(salts, ["Electrolytes/properties/salts.csv"]; verbose=verbose)

    ions = ion_groups.flattenedgroups
    _solvents = group_components(solvents)
    components = deepcopy(_solvents)
    append!(components,ions)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)

    references = String[]
    
    model = ConstW(components, _solvents, ions, isolvents, iions, value ,references)
    return model
end

function dielectric_constant(model::ConstWModel,V,T,z,_data=nothing)
    return model.value
end

is_splittable(::ConstW) = false