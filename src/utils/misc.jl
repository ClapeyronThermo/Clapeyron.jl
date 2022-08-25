"""
    vec2(x1,x2,opt=true)
Generates a correct 2-length static array `[x1,x2]`, with support for non-isbits types 
"""
function vec2(x1,x2,opt = true)
    V01,V02,_ = promote(x1,x2,opt)
    if V01 isa Base.IEEEFloat # MVector does not work on non bits types, like BigFloat
        return MVector((V01,V02))
    else
        return SizedVector{2,typeof(V01)}((V01,V02))
    end
end
"""
    dnorm(x,y,p)

Equivalent to `norm((xi-yi for (xi, yi) in zip(x, y)), p)`
"""
function dnorm(x,y,p)
    return norm((xi-yi for (xi, yi) in zip(x, y)), p)
end

include("core_utils.jl")