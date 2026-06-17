#overloads for derivative

@inline function ForwardDiff.derivative(f::WithContext{F,R1}, x::R2) where {F,R1<:Real,R2<:Real}
    T = maketagtype(f,R2)
    return ForwardDiff.extract_derivative(T, f(Dual{T}(x, one(x))))
end

@inline function ForwardDiff.derivative!(result::Union{AbstractArray,ForwardDiff.DiffResult},
                             f::WithContext{F,R1}, x::R) where {F,R1,R<:Real}
    result isa DiffResult || require_one_based_indexing(result)
    T = maketagtype(f,R)
    ydual = f(Dual{T}(x, one(x)))
    result = extract_value!(T, result, ydual)
    result = extract_derivative!(T, result, ydual)
    return result
end
