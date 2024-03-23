abstract type ZuoFirstModel <: RSPModel end

struct ZuoFirst <: ZuoFirstModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    references::Array{String,1}
end

@registermodel ZuoFirst
export ZuoFirst
function ZuoFirst(solvents,ions; userlocations::Vector{String}=String[], verbose::Bool=false)
    components = deepcopy(solvents)
    append!(components,ions)
    icomponents = 1:length(components)

    references = String[]
    
    model = ZuoFirst(components, icomponents ,references)
    return model
end

function dielectric_constant(model::ZuoFirstModel,V,T,z,_data=nothing)
    return -19.2905+29814.5/T-0.019678*T+1.318e-4*T^2-3.1144e-7*T^3
end

is_splittable(::ZuoFirst) = false