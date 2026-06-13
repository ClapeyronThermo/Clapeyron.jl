module StaticForwardDiffTagsStaticArraysExt

using StaticArrays
using ForwardDiff
using Clapeyron.StaticForwardDiffTags: maketagtype, WithContext

const ForwardDiffStatic = Base.get_extension(ForwardDiff,:ForwardDiffStaticArraysExt)

using ForwardDiff, StaticArrays
using ForwardDiff.DiffResults
using ForwardDiff: Dual,
                   gradient, hessian, jacobian, gradient!, hessian!, jacobian!,
                   extract_gradient!, extract_jacobian!, extract_value!,
                   vector_mode_gradient, vector_mode_gradient!,
                   vector_mode_jacobian, vector_mode_jacobian!, value
using DiffResults: DiffResult, ImmutableDiffResult, MutableDiffResult

const dualize = ForwardDiffStatic.dualize
const extract_gradient = ForwardDiffStatic.extract_gradient
const extract_jacobian = ForwardDiffStatic.extract_jacobian


# Gradient
@inline function ForwardDiff.vector_mode_gradient(f::F, x::StaticArray) where {F <: WithContext}
    T = maketagtype(f, eltype(x))
    return extract_gradient(T, f(dualize(T, x)), x)
end

@inline function ForwardDiff.vector_mode_gradient!(result, f::F, x::StaticArray) where {F <: WithContext}
    T = maketagtype(f, eltype(x))
    return extract_gradient!(T, result, f(dualize(T, x)))
end

# Jacobian
@inline function ForwardDiff.vector_mode_jacobian(f::F, x::StaticArray) where {F <: WithContext}
    T = maketagtype(f, eltype(x))
    return extract_jacobian(T, f(dualize(T, x)), x)
end

@inline function ForwardDiff.vector_mode_jacobian!(result, f::F, x::StaticArray) where {F <: WithContext}
    T = maketagtype(f, eltype(x))
    ydual = f(dualize(T, x))
    result = extract_jacobian!(T, result, ydual, length(x))
    result = extract_value!(T, result, ydual)
    return result
end

@inline function ForwardDiff.vector_mode_jacobian!(result::ImmutableDiffResult, f::F, x::StaticArray) where {F <: WithContext}
    T = maketagtype(f, eltype(x))
    ydual = f(dualize(T, x))
    result = DiffResults.jacobian!(result, extract_jacobian(T, ydual, x))
    result = DiffResults.value!(Base.Fix1(value, T), result, ydual)
    return result
end

# Hessian
function ForwardDiff.hessian!(result::ImmutableDiffResult, f::F, x::StaticArray) where {F <: WithContext}
    T = maketagtype(f, eltype(x))
    d1 = dualize(T, x)
    d2 = dualize(T, d1)
    fd2 = f(d2)
    val = value(T,value(T,fd2))
    grad = extract_gradient(T,value(T,fd2), x)
    hess = extract_jacobian(T,partials(T,fd2), x)
    result = DiffResults.hessian!(result, hess)
    result = DiffResults.gradient!(result, grad)
    result = DiffResults.value!(result, val)
    return result
end

end #module