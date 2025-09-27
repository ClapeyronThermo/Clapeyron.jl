#derivative logic

#ForwardDiff compiles one method per function, so separating the
#differentiation logic from the property logic allows the differentials
#to be compiled only once

"""
    ∂f∂T(model,V,T,z=SA[1.0])

Returns `∂f/∂T` at constant total volume and composition, where `f` is the total Helmholtz energy, given by `eos(model,V,T,z)`.

"""
function ∂f∂T(model,V,T,z::AbstractVector)
    f(∂T) = eos(model,V,∂T,z)
    return Solvers.derivative(f,T,∂Tag{:∂f∂T})
end

"""
    ∂f∂V(model,V,T,z)

Returns `∂f/∂V` at constant temperature `T` and composition `z`, where `f` is the total Helmholtz energy, given by `eos(model,V,T,z)`, `V` is the total volume.
"""
function ∂f∂V(model,V,T,z::AbstractVector)
    f(∂V) = a_res(model,∂V,T,z)
    ∂aᵣ∂V = Solvers.derivative(f,V,∂Tag{:∂f∂V}())
    sum(z)*Rgas(model)*T*(∂aᵣ∂V - 1/V)
end

#Returns a tuple of the form ([∂f∂V,∂f∂T],f),using the least amount of computation
"""
    ∂f(model,V,T,z)

Returns zeroth order (value) and first order derivative information of the total Helmholtz energy (given by `eos(model,V,T,z)`).
The result is given in two values:

```julia
grad_f,fval = ∂2f(model,V,T,z)
```

where:

```julia
fval   = f(V,T) = eos(model,V,T,z)

grad_f = [ ∂f/∂V; ∂f/∂T]
```

Where `V` is the total volume, `T` is the temperature and `f` is the total Helmholtz energy.
"""
function ∂f(model,V,T,z)
    f(∂V,∂T) = eos(model,∂V,∂T,z)
    _f,_df = Solvers.fgradf2(f,V,T,∂Tag{:∂f}())
    return _df,_f
end

function ∂f_vec(model,V,T,z::AbstractVector)
    _df,_f = ∂f(model,V,T,z)
    return SVector(_f,_df[1],_df[2])
end

function f∂fdV(model,V,T,z::AbstractVector)
    f(x) = eos(model,x,T,z)
    A,∂A∂V = Solvers.f∂f(f,V,∂Tag{:f∂fdV}())
    return SVector(A,∂A∂V)
end

function f∂fdT(model,V,T,z::AbstractVector)
    f(x) = eos(model,V,x,z)
    A,∂A∂T = Solvers.f∂f(f,T,∂Tag{:f∂fdT}())
    return SVector(A,∂A∂T)
end

function ∂f_res(model,V,T,z)
    f(∂V,∂T) = eos_res(model,∂V,∂T,z)
    _f,_df = Solvers.fgradf2(f,V,T,∂Tag{:f∂fdT_res}())
    return _df,_f
end

function ∂f_res_vec(model,V,T,z::AbstractVector)
    _df,_f = ∂f_res(model,V,T,z)
    return SVector(_f,_df[1],_df[2])
end

#Returns p and ∂p∂V at constant T
#it doesnt do a pass over temperature, so its
#faster that d2f when only requiring d2fdV2

"""
    p∂p∂V(model,V,T,z=SA[1.0])

Returns `p` and `∂p/∂V` at constant temperature, where `p` is the pressure = `pressure(model,V,T,z)` and `V` is the total volume.

"""
function p∂p∂V(model,V,T,z::AbstractVector=SA[1.0])
    f(∂V) = pressure(model,∂V,T,z)
    p,∂p∂V = Solvers.f∂f(f,V,∂Tag{:p∂p∂V}())
    return SVector(p,∂p∂V)
end

"""
    ∂2f(model,V,T,z)

Returns zeroth order (value), first order and second order derivative information of the total Helmholtz energy (given by `eos(model,V,T,z)`).
The result is given in three values:

```julia
hess_f,grad_f,fval = ∂2f(model,V,T,z)
```

where:
```julia
fval   = f(V,T) = eos(model,V,T,z)

grad_f = [ ∂f/∂V; ∂f/∂T]

hess_f = [ ∂²f/∂V²; ∂²f/∂V∂T
          ∂²f/∂V∂T; ∂²f/∂V²]
```

Where `V` is the total volume, `T` is the temperature and `f` is the total Helmholtz energy.
"""
function ∂2f(model,V,T,z)
    f(_V,_T) = eos(model,_V,_T,z)
    _f,_∂f,_∂2f = Solvers.∂2(f,V,T,∂Tag{:∂2f}())
    return (_∂2f,_∂f,_f)
end

"""
    ∂2p(model,V,T,z)

Returns zeroth order (value), first order and second order derivative information of the pressure.
the result is given in three values:

```julia
hess_p,grad_p,pval = ∂2p(model,V,T,z)
```

where:
```julia
pval   = p(V,T) = pressure(model,V,T,z)

grad_p = [ ∂p/∂V; ∂p/∂T]

hess_p = [ ∂²p/∂V²; ∂²p/∂V∂T
          ∂²p/∂V∂T; ∂²p/∂V²]
```

Where `V` is the total volume, `T` is the temperature and `p` is the pressure.
"""
function ∂2p(model,V,T,z)
    f(_V,_T) = pressure(model,_V,_T,z)
    _f,_∂f,_∂2f = Solvers.∂2(f,V,T,∂Tag{:∂2p}())
    return (_∂2f,_∂f,_f)
end

"""
    f_hess(model,V,T,z)

Returns the second order volume `V` and temperature `T` derivatives of the total Helmholtz energy `f` (given by `eos(model,V,T,z)`). The result is given in a 2x2 `SMatrix`, in the form:

```julia
[ ∂²f/∂V²  ∂²f/∂V∂T
 ∂²f/∂V∂T  ∂²f/∂T²]
```

Use this instead of the `∂2f` if you only need second order information. `∂2f` also gives zeroth and first order derivative information, but due to a bug in the used AD, it allocates more than necessary.
"""
function f_hess(model,V,T,z)
    f(w) = eos(model,first(w),last(w),z)
    V,T = promote(V,T)
    VT_vec = SVector(V,T)
    return Solvers.hessian(f,VT_vec)
end

"""
    p∂p∂2p(model,V,T,z=SA[1.0])

Returns the pressure `p` and their first and second volume derivatives `∂p/∂V` and `∂²p/∂V²`, in a single ForwardDiff pass.

"""
function p∂p∂2p(model,V,T,z=SA[1.0])
    f(∂V) = pressure(model,∂V,T,z)
    _, ∂²A∂V², ∂³A∂V³ = Solvers.f∂f∂2f(f,V,∂Tag{:p∂p∂2p}())
    return ∂²A∂V², ∂³A∂V³
end

"""
    ∂²f∂T²(model,V,T,z=SA[1.0])

Returns `∂²A/∂T²` via Autodiff. Used mainly for ideal gas properties. It is recommended to overload this function for ideal models, as is equivalent to -Cv(T)/T.

"""
function ∂²f∂T²(model,V,T,z)
    A(x) = eos(model,V,x,z)
    _,_,∂²A∂T² = Solvers.f∂f∂2f(A,V,∂Tag{:∂²f∂T²}())
    return ∂²A∂T²
end

#derivarive logic: model Dual numbers:

#as of Clapeyron 0.6.10, there is limited support for using models with dual numbers
#PCSAFT, sPCSAFT, SAFTVRMie, SAFTVRMie15 support using dual numbers, (and any other number type)
#for iterative methods, it is more efficient to reconstruct the model with the primal value instead of the full value.

function Solvers.primalval(model::EoSModel)
    return _primalval(model,eltype(model))
end

function _primalval(model::EoSModel,::Type{T}) where T <: ForwardDiff.Dual
    return Solvers.primalval_struct(model)
end

_primalval(model::EoSModel,::T) where T = model