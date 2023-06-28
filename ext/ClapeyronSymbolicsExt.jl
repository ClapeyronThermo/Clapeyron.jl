module ClapeyronSymbolicsExt
using  Clapeyron
using Symbolics

using Clapeyron:log1p,log,sqrt

Symbolics.@register_symbolic Clapeyron.Solvers.log(x)
Symbolics.@register_symbolic Clapeyron.Solvers.log1p(x)
Symbolics.@register_symbolic Clapeyron.Solvers.sqrt(x)

Symbolics.derivative(::typeof(log), args::NTuple{1,Any}, ::Val) where {N} = 1/args[1]
Symbolics.derivative(::typeof(log1p), args::NTuple{1,Any}, ::Val) where {N} = 1/(1 + args[1])
Symbolics.derivative(::typeof(sqrt), args::NTuple{1,Any}, ::Val) where {N} = 1 / (2sqrt(args[1]))

function Solvers.derivative(f::F,x::Symbolics.Num) where {F}
    fx = f(x)
    dfx = Symbolics.derivative(fx,x)
    Symbolics.simplify(dfx)
end

function Solvers.gradient(f::F,x::A) where {F,A<:AbstractArray{N},N::Symbolics.Num}
    fx = f(x)
    dfx = Symbolics.gradient(fx,x)
    Symbolics.simplify(dfx)
end

function Solvers.hessian(f::F,x::A) where {F,A<:AbstractArray{N},N::Symbolics.Num}
    fx = f(x)
    dfx = Symbolics.gradient(fx,x)
    Symbolics.simplify(dfx)
end

end #module