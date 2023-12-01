abstract type ConstRSPModel <: RSPModel end

struct ConstRSP <: ConstRSPModel
    components::Array{String,1}
    neutral::Array{String,1}
    ions::Array{String,1}
    icomponents::UnitRange{Int}
    ineutral::UnitRange{Int}
    iions::UnitRange{Int}
    value::Float64
    references::Array{String,1}
end

@registermodel ConstRSP
export ConstRSP
function ConstRSP(solvents,ions; userlocations::Vector{String}=String[], value =  78.4, verbose::Bool=false)
    components = deepcopy(solvents)
    append!(components,ions)
    isolvents = 1:length(solvents)
    iions = (length(solvents)+1):length(components)
    icomponents = 1:length(components)

    references = String[]
    
    model = ConstRSP(components, solvents, ions, icomponents, isolvents, iions, value ,references)
    return model
end

function dielectric_constant(model::ConstRSPModel,V,T,z,_data=nothing)
    return model.value
end

is_splittable(::ConstRSP) = false