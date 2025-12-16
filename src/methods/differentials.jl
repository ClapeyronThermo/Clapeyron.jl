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
    return Solvers.derivative(f,T)
end

"""
    ∂f∂V(model,V,T,z)

Returns `∂f/∂V` at constant temperature `T` and composition `z`, where `f` is the total Helmholtz energy, given by `eos(model,V,T,z)`, `V` is the total volume.
"""
function ∂f∂V(model,V,T,z::AbstractVector)
    f(∂V) = a_res(model,∂V,T,z)
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

grad_f = [ ∂f/∂V; ∂f/∂T]
```

Where `V` is the total volume, `T` is the temperature and `f` is the total Helmholtz energy.
"""
function ∂f(model,V,T,z)
    f(∂V,∂T) = eos(model,∂V,∂T,z)
    _f,_df = Solvers.fgradf2(f,V,T)
    return _df,_f
end

function ∂f_vec(model,V,T,z::AbstractVector)
    _df,_f = ∂f(model,V,T,z)
    return SVector(_f,_df[1],_df[2])
end

function f∂fdV(model,V,T,z::AbstractVector)
    f(x) = eos(model,x,T,z)
    A,∂A∂V = Solvers.f∂f(f,V)
    return SVector(A,∂A∂V)
end

function f∂fdT(model,V,T,z::AbstractVector)
    f(x) = eos(model,V,x,z)
    A,∂A∂T = Solvers.f∂f(f,T,)
    return SVector(A,∂A∂T)
end

function ∂f_res(model,V,T,z)
    f(∂V,∂T) = eos_res(model,∂V,∂T,z)
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

"""
    p∂p∂V(model,V,T,z=SA[1.0])

Returns `p` and `∂p/∂V` at constant temperature, where `p` is the pressure = `pressure(model,V,T,z)` and `V` is the total volume.

"""
function p∂p∂V(model,V,T,z::AbstractVector=SA[1.0])
    f(∂V) = pressure(model,∂V,T,z)
    p,∂p∂V = Solvers.f∂f(f,V)
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
    _f,_∂f,_∂2f = Solvers.∂2(f,V,T)
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
    p, ∂²A∂V², ∂³A∂V³ = Solvers.f∂f∂2f(f,V)
    return p, ∂²A∂V², ∂³A∂V³
end

"""
    ∂²f∂T²(model,V,T,z=SA[1.0])

Returns `∂²A/∂T²` via Autodiff. Used mainly for ideal gas properties. It is recommended to overload this function for ideal models, as is equivalent to -Cv(T)/T.

"""
function ∂²f∂T²(model,V,T,z)
    A(_T) = eos(model,V,_T,z)
    _,_,∂²A∂T² = Solvers.f∂f∂2f(A,T)
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

"""
    __gradients_for_root_finders(x::AbstractVector{T},tups::Tuple,tups_primal::Tuple,f::Function) where T<:Real

Computes the gradients of `x` with respect to the relevant parameters in `tups` under the condition that `x` is implicitly defined through the root finding problem `f(x,tups) = 0`. 
The function uses the implicit function theorem to compute the gradients efficiently through the reconstruction of Duals.

!!! note 
    Currently only supports first order AD. Trying to differentiate nested Duals or Duals with different tags will throw an error.
"""
function __gradients_for_root_finders(x::AbstractVector{T},tups::T1,tups_primal::T2,f::F) where {T<:Real,T1,T2,F}
    # check and return primal of no duals
    !any(has_dual,tups) && return x

    # Checks
    implicit_ad_check(tups)

    # compute partials
    ∂f_∂θ_dual = f(x,tups)

    #guard against NaN in input or output
    if any(isnan,x) || any(isnan,∂f_∂θ_dual)
        _0 = zero(eltype(∂f_∂θ_dual))
        nan = _0/_0
        return nan .* x
    end

    return __gradients_for_root_finders_solver(f,x,tups_primal,∂f_∂θ_dual)
end

function __gradients_for_root_finders(x::T,tups::T1,primal_tups::T2,f::F) where {T<:Real,T1,T2,F}
    x_vec = SVector(x)
    f_vec(x_vec,tupsx) = SVector(f(x_vec[1],tupsx))
    x_dual_vec = __gradients_for_root_finders(x_vec,tups,primal_tups,f_vec)
    return only(x_dual_vec)
end

__gradients_for_root_finders(::Union{AbstractArray{T},T},_,_,_) where T<:ForwardDiff.Dual = error("Input `x` cannot be a dual")

function __partial_to_svec(x::ForwardDiff.Dual{T,V,N}) where {T,V,N}
    return SVector{N,V}(ForwardDiff.partials(x))
end

function __gradients_for_root_finders_solver(f::F,x::AbstractVector{V1},tups_primal,∂f_∂θ_dual::AbstractVector{ForwardDiff.Dual{TAG,V2,Npartials}}) where {F,V1,TAG,V2,Npartials}

    ∂f_∂x = DifferentiationInterface.jacobian(f,DifferentiationInterface.AutoForwardDiff(),x,DifferentiationInterface.Constant(tups_primal))
    LinearAlgebra.checksquare(∂f_∂x)
    # get dual
    ∂f_∂θi = map(__partial_to_svec,∂f_∂θ_dual)
    ∂f_∂θ = transpose(reduce(hcat,∂f_∂θi))

    # gradient through implicit function theorem
    dx_dθ = - ∂f_∂x \ ∂f_∂θ

    to_dual(x1,x2) = ForwardDiff.Dual{TAG,V2,Npartials}(x1,ForwardDiff.Partials(NTuple{Npartials,V2}(x2)))
    x_dual2 = map(to_dual,x,eachrow(dx_dθ))
    return x_dual2
end

function nested_ad_check(a::A) where A
    AT = eltype(a)
    if AT <: ForwardDiff.Dual
        V = ForwardDiff.valtype(AT)
        V isa ForwardDiff.Dual && throw(NestedADError("Found nested Duals of type $AT. This is currently not supported in implicit differentiation."))
    end
    return nothing
end

nested_ad_check(a::Tuple) = foreach(nested_ad_check,a)

function multiple_tag_ad_check(a::T) where T <: Tuple
    f(x) = eltype(x) <: ForwardDiff.Dual
    n_dual = count(f,a)
    if n_dual == 1
        return nothing
    end

    #TODO write as @generated so it does not allocate
    tag(x) = ForwardDiff.tagtype(eltype(x))
    duals = findall(map(f,a))
    tags = map(tag,a)
    valid_tags = tags[duals]
    t1 = first(valid_tags)
    t = Base.tail(valid_tags)
    for ti in t
        if t1 != ti
            msg = "Found multiple Dual tags: $t1 and $ti. This is currently not supported in implicit differentiation."
            throw(MultipleTagError(msg))
        end
    end
    return nothing
end

multiple_tag_ad_check(a) = nothing

function implicit_ad_check(a)
    nested_ad_check(a)
end

function implicit_ad_check(a::Tuple)
    nested_ad_check(a)
    multiple_tag_ad_check(a)
end