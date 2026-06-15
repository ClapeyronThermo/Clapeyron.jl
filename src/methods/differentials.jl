#derivative logic

#ForwardDiff compiles one method per function, so separating the
#differentiation logic from the property logic allows the differentials
#to be compiled only once

struct ∂₁f end

"""
    ∂f∂T(model,V,T,z=SA[1.0])

Returns `∂f/∂T` at constant total volume `V` and composition `z`, where `f` is the total Helmholtz energy, given by `eos(model,V,T,z)`.

"""
function ∂f∂T(model,V,T,z::AbstractVector)
    f = @deferred_T(eos,∂₁f)
    return Solvers.derivative(f,T)
end

"""
    ∂f∂V(model,V,T,z)

Returns `∂f/∂V` at constant temperature `T` and composition `z`, where `f` is the total Helmholtz energy, given by `eos(model,V,T,z)`, `V` is the total volume.
"""
function ∂f∂V(model,V,T,z::AbstractVector)
    f = @deferred_V(a_res,∂₁f)
    ∂aᵣ∂V = Solvers.derivative(f,V)
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

grad_f = [∂f/∂V; ∂f/∂T]
```

Where `V` is the total volume, `T` is the temperature and `f` is the total Helmholtz energy.
"""
function ∂f(model,V,T,z)
    f = @deferred_VT(eos,∂₁f)
    _f,_df = Solvers.fgradf2(f,V,T)
    return _df,_f
end

function ∂f_vec(model,V,T,z::AbstractVector)
    _df,_f = ∂f(model,V,T,z)
    return SVector(_f,_df[1],_df[2])
end

function f∂fdV(model,V,T,z::AbstractVector)
    f = @deferred_V(eos,∂₁f)
    A,∂A∂V = Solvers.f∂f(f,V)
    return SVector(A,∂A∂V)
end

function f∂fdV_res(model,V,T,z::AbstractVector)
    f = @deferred_V(eos_res,∂₁f)
    Ar,∂Ar∂V = Solvers.f∂f(f,V)
    return SVector(Ar,∂Ar∂V)
end

function f∂fdT(model,V,T,z::AbstractVector)
    f = @deferred_T(eos,∂₁f)
    A,∂A∂T = Solvers.f∂f(f,T)
    return SVector(A,∂A∂T)
end

function ∂f_res(model,V,T,z)
    f = @deferred_VT(eos_res,∂₁f)
    _f,_df = Solvers.fgradf2(f,V,T)
    return _df,_f
end

function ∂f_res_vec(model,V,T,z::AbstractVector)
    _df,_f = ∂f_res(model,V,T,z)
    return SVector(_f,_df[1],_df[2])
end

#Returns p and ∂p∂V at constant T
#it doesnt do a pass over temperature, so its
#faster that d2f when only requiring d2fdV2

struct ∂₁p end

"""
    p∂p∂V(model,V,T,z=SA[1.0])

Returns `p` and `∂p/∂V` at constant temperature `T`, where `p` is the pressure = `pressure(model,V,T,z)` and `V` is the total volume.

"""
function p∂p∂V(model,V,T,z::AbstractVector=SA[1.0])
    f = @deferred_V(pressure,∂₁p)
    p,∂p∂V = Solvers.f∂f(f,V)
    return SVector(p,∂p∂V)
end

function p∂p∂rho(model, rho, T, z=SA[1.0])
    n   = sum(z)
    V   = n / rho
    p, dpdV = p∂p∂V(model, V, T, z)
    dVdρ    = -V / rho
    d2Vdρ2  =  2 * V / (rho*rho)
    dpdrho   = dpdV * dVdρ
    return SVector(p, dpdrho)
end

function ∂p∂rho(model, rho, T, z=SA[1.0])
    _,dpdrho = p∂p∂rho(model,rho,T,z)
    return dpdrho
end


"""
    ∂p∂T(model,V,T,z=SA[1.0])

Returns `∂p/∂T` at constant temperature `T`, where `p` is the pressure = `pressure(model,V,T,z)` and `V` is the total volume.

"""
function ∂p∂T(model,V,T,z::AbstractVector=SA[1.0])
    f = @deferred_T(pressure,∂₁p)
    return Solvers.derivative(f,T)
end

struct ∂₂f end
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
    f = @deferred_VT(eos,∂₂f)
    _f,_∂f,_∂2f = Solvers.∂2(f,V,T)
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
    f = @deferred_VT(eos,∂₂f)
    V,T = promote(V,T)
    VT_vec = SVector(V,T)
    return Solvers.hessian(f,VT_vec)
end

"""
    ∂²f∂T²(model,V,T,z=SA[1.0])

Returns `∂²A/∂T²` via Autodiff. Used mainly for ideal gas properties. It is recommended to overload this function for ideal models, as is equivalent to -Cv(T)/T.

"""
function ∂²f∂T²(model,V,T,z)
    f = @deferred_T(eos,∂₂f)
    _,_,∂²A∂T² = Solvers.f∂f∂2f(f,T)
    return ∂²A∂T²
end

struct ∂₂p end

"""
    ∂2p(model,V,T,z)

Returns zeroth order (value), first order and second order derivative information of the pressure.
The result is given in three values:

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
    f = @deferred_VT(pressure,∂₂p)
    _f,_∂f,_∂2f = Solvers.∂2(f,V,T)
    return (_∂2f,_∂f,_f)
end

function ∂2p_vec(model,V,T,z)
    _∂2f,_∂f,_f = ∂2p(model,V,T,z)
    return SVector(_f,_∂f[1],_∂f[2],_∂2f[1,1],_∂2f[2,2],_∂2f[1,2])
end

"""
    p∂p∂2p(model,V,T,z=SA[1.0])

Returns the pressure `p` and their first and second volume derivatives `∂p/∂V` and `∂²p/∂V²`, in a single ForwardDiff pass.

"""
function p∂p∂2p(model,V,T,z=SA[1.0])
    f = @deferred_V(pressure,∂₂p)
    p, ∂²A∂V², ∂³A∂V³ = Solvers.f∂f∂2f(f,V)
    return SVector(p, ∂²A∂V², ∂³A∂V³)
end

function p∂p∂2p_rho(model, rho, T, z=SA[1.0])
    n   = sum(z)
    V   = n / rho
    #rho = n / V, 1/rho = V/n
    p, dpdV, d2pdV2 = p∂p∂2p(model, V, T, z)

    # ρ = n/V  →  V = n/ρ
    # dV/dρ   = -V²/n
    # d²V/dρ² =  2V³/n²
    dVdρ    = -V / rho
    d2Vdρ2  =  2 * V / (rho*rho)
    # Chain rule
    dpdrho   = dpdV * dVdρ
    d2pdrho2 = d2pdV2 * dVdρ*dVdρ + dpdV * d2Vdρ2
    return SVector(p, dpdrho, d2pdrho2)
end

#derivative logic: model Dual numbers:

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

# use IFTDuals: ift for implicit differentiation
"""
    __gradients_for_root_finders(x::AbstractVector{T},tups::Tuple,tups_primal::Tuple,f::Function) where T<:Real

Computes the gradients of `x` with respect to the relevant parameters in `tups` under the condition that `x` is implicitly defined through the root finding problem `f(x,tups) = 0`.
The function uses the implicit function theorem to compute the gradients efficiently through the reconstruction of Duals. We use the IFTDuals.jl package for this purpose, which has some restrictions, currently mixed nested Duals (i.e. different tags) are not supported.
"""
function __gradients_for_root_finders(x::Union{AbstractArray{T},T},tups,tups_primal,f::Function) where T<:Real # tups not restructed to Tuple
    if any(isnan,x) # guard against NaN in input, do not need Dual types here?
        return x
    end
    return ift(x,f,tups,tups_primal) # use IFTDuals package, returns primal if tups has no duals
end

IFTDuals.promote_my_type(m::EoSModel) = eltype(m)
