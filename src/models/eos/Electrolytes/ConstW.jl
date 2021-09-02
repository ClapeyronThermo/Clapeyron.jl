abstract type RSPModel <: EoSModel end
abstract type ConstWModel <: RSPModel end

struct ConstWParam <: EoSParam
end

@newmodelsimple ConstW ConstWModel ConstWParam

export ConstW
function ConstW(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    model = ConstW(ConstWParam())
    return model
end

function RSP(model::ConstWModel,V,T,z,)
    return 78.4
end

is_splittable(::ConstW) = false