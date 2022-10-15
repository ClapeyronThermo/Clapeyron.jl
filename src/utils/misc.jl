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

#this is never used in a critical path, so we just use a default copying method
if VERSION >= v"1.7"
    keepat!(a,inds) = Base.keepat!(a,inds)
else
    function keepat!(a,inds)
        b = a[inds]
        resize!(a,length(b))
        a .= b
        return a
    end
end

include("core_utils.jl")