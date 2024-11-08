#derivative logic

#ForwardDiff compiles one method per function, so separating the
#differentiation logic from the property logic allows the differentials
#to be compiled only once

"""
    ∂f∂T(model,V,T,z=SA[1.0])

returns `∂f/∂T` at constant total volume and composition, where f is the total helmholtz energy, given by `eos(model,V,T,z)`

"""
function ∂f∂T(model,V,T,z)
    f(∂T) = eos(model,V,∂T,z)
    return Solvers.derivative(f,T)
end

"""
    ∂f∂V(model,V,T,z)

returns `∂f/∂V` at constant temperature and composition, where f is the total helmholtz energy, given by `eos(model,V,T,z)`, and V is the total volume
"""
function ∂f∂V(model,V,T,z)
    f(∂V) = a_res(model,∂V,T,z)
    ∂aᵣ∂V = Solvers.derivative(f,V)
    sum(z)*Rgas(model)*T*(∂aᵣ∂V - 1/V)
end

#returns a tuple of the form ([∂f∂V,∂f∂T],f),using the least amount of computation
"""
    ∂f(model,V,T,z)

returns zeroth order (value) and first order derivative information of the total helmholtz energy (given by `eos(model,V,T,z)`).
the result is given in two values:

```julia
grad_f,fval = ∂2f(model,V,T,z)
```

where:

```julia
fval   = f(V,T) = eos(model,V,T,z)

grad_f = [ ∂f/∂V; ∂f/∂T]

```

Where `V` is the total volume, `T` is the temperature and `f` is the total helmholtz energy.
"""
function ∂f(model,V,T,z)
    f(∂V,∂T) = eos(model,∂V,∂T,z)
    _f,_df = Solvers.fgradf2(f,V,T)
    return _df,_f
end

function ∂f_vec(model,V,T,z)
    _df,_f = ∂f(model,V,T,z)
    return SVector(_f,_df[1],_df[2])
end

function f∂fdV(model,V,T,z)
    f(x) = eos(model,x,T,z)
    A,∂A∂V = Solvers.f∂f(f,V)
    return SVector(A,∂A∂V)
end

function f∂fdT(model,V,T,z)
    f(x) = eos(model,V,x,z)
    A,∂A∂T = Solvers.f∂f(f,T)
    return SVector(A,∂A∂T)
end

function ∂f_res(model,V,T,z)
    f(∂V,∂T) = eos_res(model,∂V,∂T,z)
    _f,_df = Solvers.fgradf2(f,V,T)
    return _df,_f
end
#returns p and ∂p∂V at constant T
#it doesnt do a pass over temperature, so its
#faster that d2f when only requiring d2fdV2

"""
    p∂p∂V(model,V,T,z=SA[1.0])

returns `p` and `∂p/∂V` at constant temperature, where p is the pressure = `pressure(model,V,T,z)` and `V` is the total Volume.

"""
function p∂p∂V(model,V,T,z=SA[1.0])
    f(∂V) = pressure(model,∂V,T,z)
    p,∂p∂V = Solvers.f∂f(f,V)
    return SVector(p,∂p∂V)
end

"""
    ∂2f(model,V,T,z)

returns zeroth order (value), first order and second order derivative information of the total helmholtz energy (given by `eos(model,V,T,z)`).
the result is given in three values:

```
hess_f,grad_f,fval = ∂2f(model,V,T,z)
```

where:
```
fval   = f(V,T) = eos(model,V,T,z)

grad_f = [ ∂f/∂V; ∂f/∂T]

hess_f = [ ∂²f/∂V²; ∂²f/∂V∂T
          ∂²f/∂V∂T; ∂²f/∂V²]
 ```

Where `V` is the total volume, `T` is the temperature and `f` is the total helmholtz energy.
"""
function ∂2f(model,V,T,z)
    f(_V,_T) = eos(model,_V,_T,z)
    _f,_∂f,_∂2f = Solvers.∂2(f,V,T)
    return (_∂2f,_∂f,_f)
end

"""
    ∂2p(model,V,T,z)

returns zeroth order (value), first order and second order derivative information of the pressure.
the result is given in three values:

```
hess_p,grad_p,pval = ∂2p(model,V,T,z)
```

where:
```
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

returns the second order volume (`V`) and temperature (`T`) derivatives of the total helmholtz energy (given by `eos(model,V,T,z)`). the result is given in a 2x2 `SMatrix`, in the form:

```
[ ∂²f/∂V²  ∂²f/∂V∂T
 ∂²f/∂V∂T  ∂²f/∂T²]
 ```

use this instead of the ∂2f if you only need second order information. ∂2f also gives zeroth and first order derivative information, but due to a bug in the used AD, it allocates more than necessary.
"""
function f_hess(model,V,T,z)
    f(w) = eos(model,first(w),last(w),z)
    V,T = promote(V,T)
    VT_vec = SVector(V,T)
    return Solvers.hessian(f,VT_vec)
end

"""
    ∂²³f(model,V,T,z=SA[1.0])

returns `∂²A/∂V²` and `∂³A/∂V³`, in a single ForwardDiff pass. used mainly in `crit_pure` objective function

"""
function ∂²³f(model,V,T,z=SA[1.0])
    f(∂V) = pressure(model,∂V,T,z)
    _, ∂²A∂V², ∂³A∂V³ = Solvers.f∂f∂2f(f,V)
    return ∂²A∂V², ∂³A∂V³
end

"""
    ∂²f∂T²(model,V,T,z=SA[1.0])

returns `∂²A/∂T²` via Autodiff. Used mainly for ideal gas properties. It is recommended to overload this function for ideal models, as is equivalent to -Cv(T)/T

"""
function ∂²f∂T²(model,V,T,z)
    A(x) = eos(model,V,x,z)
    ∂A∂T(x) = Solvers.derivative(A,x)
    ∂²A∂T²(x) = Solvers.derivative(∂A∂T,x)
    return ∂²A∂T²(T)
end

function d2fdt2(model,V,T,z)
    A(x) = eos(model,V,x,z)
    ∂A∂T(x) = Solvers.derivative(A,x)
    ∂²A∂T²(x) = Solvers.derivative(∂A∂T,x)
    return ∂²A∂T²(T)
end


const _d23f = ∂²³f