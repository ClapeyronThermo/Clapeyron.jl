module ClapeyronSymbolicsExt
using Clapeyron
using Clapeyron.ForwardDiff
using Clapeyron.Solvers
using Clapeyron.EoSFunctions

using Symbolics

using Clapeyron: log1p,log,sqrt,^
using Clapeyron: SA

Solvers.log(x::Num) = Base.log(x)
Solvers.log1p(x::Num) = Base.log1p(x)
Solvers.sqrt(x::Num) = Base.log1p(x)
EoSFunctions.xlogx(x::Num,k) = x*Base.log(x*k)
EoSFunctions.xlogx(x::Num) = x*Base.log(x)
Solvers.:^(x::Num,y) = Base.:^(x,y)
Solvers.:^(x,y::Num) = Base.:^(x,y)
Solvers.:^(x::Num,y::Int) = Base.:^(x,y)
Solvers.:^(x::Num,y::Num) = Base.:^(x,y)

function Solvers.derivative(f::F,x::Symbolics.Num) where {F}
    fx = f(x)
    dfx = Symbolics.derivative(fx,x)
    Symbolics.simplify(dfx)
end

function Solvers.gradient(f::F,x::A) where {F,A<:AbstractArray{Symbolics.Num}}
    fx = f(x)
    Symbolics.gradient(fx,x)
end

function Solvers.hessian(f::F,x::A) where {F,A<:AbstractArray{Symbolics.Num}}
    fx = f(x)
    Symbolics.hessian(fx,x)
end
#=
function Solvers.:^(x::Num, y::ForwardDiff.Dual{Ty}) where Ty
    _y = ForwardDiff.value(y)
    fy = x^_y
    dy = log(x)*fy
    ForwardDiff.Dual{Ty}(fy, dy * ForwardDiff.partials(y))
end

function Solvers.:^(x::ForwardDiff.Dual{Tx},y::Symbolics.Num) where Tx
    _x = ForwardDiff.value(x)
    fx = _x^y
    dx = y*_x^(y-1)
    #partials(x) * y * ($f)(v, y - 1)
    ForwardDiff.Dual{Tx}(fx, dx * ForwardDiff.partials(x))
end

function Solvers.:^(x::ForwardDiff.Dual{Tx,Num},y::Int) where Tx
    _x = ForwardDiff.value(x)
    fx = _x^y
    dx = y*_x^(y-1)
    #partials(x) * y * ($f)(v, y - 1)
    ForwardDiff.Dual{Tx}(fx, dx * ForwardDiff.partials(x))
end
=#

function Solvers.f∂f(f::F, x::Num) where {F}
    fx = f(x)
    return fx,Symbolics.derivative(fx,x)
end

function Solvers.f∂f∂2f(f::F, x::Num) where {F}
    fx,dfx = Solvers.f∂f(f,x)
    return fx,dfx,Symbolics.derivative(dfx,x)
end

function fgradf2_sym(f,V,T)
    fvt = f(V,T)
    dv = Symbolics.derivative(fvt,V)
    dT = Symbolics.derivative(fvt,T)
    return fvt,SA[dv,dT]
end

function Solvers.fgradf2(f,V::Num,T)
    @variables T̃
    fvt,dvt = fgradf2_sym(f,V,T̃)
    t_dict = Dict(T̃ => T)
    fvt,Symbolics.substitute(dvt,t_dict)
end

function Solvers.fgradf2(f,V,T::Num)
    @variables Ṽ
    fvt,dvt = fgradf2_sym(f,Ṽ,T)
    v_dict = Dict(Ṽ => V)
    fvt,Symbolics.substitute(dvt,v_dict)
end

Solvers.fgradf2(f,V::Num,T::Num) = fgradf2_sym(f,V,T)

gradient2_sym(f,V,T) = last(fgradf2_sym(f,V,T))

function Solvers.gradient2(f,V::Num,T)
    @variables T̃
    dvt = gradient2_sym(f,V,T̃)
    t_dict = Dict(T̃ => T)
    Symbolics.substitute(dvt,t_dict)
end

function Solvers.gradient2(f,V,T::Num)
    @variables Ṽ
    dvt = gradient2_sym(f,Ṽ,T)
    v_dict = Dict(Ṽ => V)
    Symbolics.substitute(dvt,v_dict)
end

Solvers.gradient2(f,V::Num,T::Num) = gradient2_sym(f,V,T)

function ∂2_sym(f,V,T)
    fvt,gvt = fgradf2_sym(f,V,T)
    x = SA[V,T]
    hvt = Symbolics.jacobian(gvt,x)
    return fvt,gvt,hvt
end

function Solvers.∂2(f,V::Num,T)
    @variables T̃
    fvt,gvt,hvt = ∂2_sym(f,V,T̃)
    t_dict = Dict(T̃ => T)
    fvt,Symbolics.substitute(gvt,t_dict),Symbolics.substitute(hvt,t_dict)
end

function Solvers.∂2(f,V,T::Num)
    @variables Ṽ
    fvt,gvt,hvt = ∂2_sym(f,Ṽ,T)
    v_dict = Dict(Ṽ => V)
    fvt,Symbolics.substitute(gvt,v_dict),Symbolics.substitute(hvt,v_dict)
end

Solvers.∂2(f,V::Num,T::Num) = ∂2_sym(f,V,T)

SymbolicUtils.promote_symtype(::typeof(Clapeyron.volume), args...) = Real
SymbolicUtils.promote_symtype(::typeof(Clapeyron.eos), args...) = Real 
#SymbolicUtils.promote_symtype(::typeof(Clapeyron.dfdv), args...) = Real 

const SymReal = Union{Num,SymbolicUtils.Symbolic{<:Real}}
const SymZ = Union{AbstractVector{Num},Symbolics.SymArray,Num,SymbolicUtils.Symbolic{<:Real},Symbolics.Arr{Num, 1}}
const RealZ = Union{AbstractVector{<:Real},Real}

@register_symbolic Clapeyron._volume(model::EoSModel,p,T,arr::AbstractVector,
                                     sym::Union{String,Symbol},bool::Bool,x::Union{Real,Nothing})
@register_symbolic Clapeyron._volume(model::EoSModel,p,T,arr::Symbolics.Arr{Num,1},
                                     sym::Union{String,Symbol},bool::Bool,x::Union{Real,Nothing})
#@register_symbolic Clapeyron._volume(model::EoSModel,p,T,z::Symbolics.Arr{Num,1},phase::Symbol,threaded::Bool,vol0::Union{Real,Nothing}) false
#@register_symbolic Clapeyron._volume(model::EoSModel,p,T,z::AbstractVector,phase::Symbol,threaded::Bool,vol0::Union{Real,Nothing}) false
@register_symbolic Clapeyron.eos(model::EoSModel,p,T,z::AbstractVector) false
@register_symbolic Clapeyron.eos(model::EoSModel,p,T,z::Symbolics.Arr{Num,1}) false
@register_symbolic Clapeyron.eos(model::EoSModel,p,T) false
Symbolics.@register_array_symbolic Clapeyron.∂f_vec(model::EoSModel,p,T,z::AbstractVector) begin
size=(3,)
end

Symbolics.@register_array_symbolic Clapeyron.∂f_vec(model::EoSModel,p,T,z::Symbolics.Arr{Num,1}) begin
size=(3,)
end

@register_symbolic Clapeyron.∂f∂V(model::EoSModel,V,T,z::Symbolics.Arr{Num,1})
@register_symbolic Clapeyron.∂f∂V(model::EoSModel,V,T,z::AbstractVector) false
@register_symbolic Clapeyron.∂f∂V(model::Clapeyron.SecondVirialModel,V,T,z::Symbolics.Arr{Num,1}) false
@register_symbolic Clapeyron.∂f∂V(model::Clapeyron.SecondVirialModel,V,T,z::AbstractVector) false
@register_symbolic Clapeyron.∂f∂T(model::EoSModel,V,T,z::Symbolics.Arr{Num,1})
@register_symbolic Clapeyron.∂f∂T(model::EoSModel,V,T,z::AbstractVector) false

Symbolics.@register_array_symbolic Clapeyron.f∂fdV(model::EoSModel,p,T,z::AbstractVector) begin
    size=(2,)
end
Symbolics.@register_array_symbolic Clapeyron.f∂fdV(model::EoSModel,p,T,z::Symbolics.Arr{Num,1}) begin
    size=(2,)
end

Symbolics.@register_array_symbolic Clapeyron.VT_partial_property(model::EoSModel,V,T,z::AbstractVector,property::Function) begin
    size = size(z)
    eltype = Base.promote_eltype(model,V,T,z)
end

Symbolics.@register_array_symbolic Clapeyron.VT_partial_property(model::EoSModel,V,T,z::Symbolics.Arr{Num,1},property::Function) begin
    size = size(z)
    eltype = Base.promote_eltype(model,V,T,z)
end

Symbolics.@register_array_symbolic Clapeyron.VT_partial_property(model::EoSModel,V,T,z::AbstractVector,property::typeof(Clapeyron.volume)) begin
    size = size(z)
    eltype = Base.promote_eltype(model,V,T,z)
end

Symbolics.@register_array_symbolic Clapeyron.VT_partial_property(model::EoSModel,V,T,z::Symbolics.Arr{Num,1},property::typeof(Clapeyron.volume)) begin
    size = size(z)
    eltype = Base.promote_eltype(model,V,T,z)
end

Symbolics.@register_array_symbolic Clapeyron.p∂p∂V(model::EoSModel,p,T,z::AbstractVector) begin
    size=(2,)
end
Symbolics.@register_array_symbolic Clapeyron.p∂p∂V(model::EoSModel,p,T,z::Symbolics.Arr{Num,1}) begin
    size=(2,)
end

end #module