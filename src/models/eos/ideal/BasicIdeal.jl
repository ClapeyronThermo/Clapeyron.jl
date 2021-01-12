struct BasicIdealParam <: EoSParam
end

abstract type BasicIdealModel <: IdealModel end
struct BasicIdeal <: BasicIdealModel
    modelname::String
    params::BasicIdealParam
    BasicIdeal(params::BasicIdealParam) = new("BasicIdeal", params)
end

function BasicIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    return BasicIdeal(BasicIdealParam())
end

function BasicIdeal(; userlocations::Array{String,1}=String[], verbose=false)
    return BasicIdeal(BasicIdealParam())
end

function a_ideal(model::BasicIdeal, V, T, z)
    x = z/∑(z)
    return ∑(x .* log.(z/V)) - 1
end
