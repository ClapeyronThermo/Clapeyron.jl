struct BasicIdealParam <: EoSParam
end

abstract type BasicIdealModel <: IdealModel end
struct BasicIdeal <: BasicIdealModel
    params::BasicIdealParam
end


export BasicIdeal
function BasicIdeal(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false)
    return BasicIdeal(BasicIdealParam())
end

function BasicIdeal(; userlocations::Array{String,1}=String[], verbose=false)
    return BasicIdeal(BasicIdealParam())
end
is_splittable(::BasicIdeal) = false

function a_ideal(model::BasicIdeal, V, T, z)
    N = ∑(z)
    #x = z/∑(z)
    res = ∑(xlogx,z) 
    res /= N 
    res -= log(V) 
    res -= 1.5*log(T)
    res -= one(res)
    # ∑(x .* log.(z/V)) - 1 original formulation, prone no NaN when passing pure Fractions
    return res
end
