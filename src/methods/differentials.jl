#derivative logic

#ForwardDiff compiles one method per function, so separating the
#differentiation logic from the property logic allows the differentials
#to be compiled only once


"""
    ∂f∂T(model,V,T,z=SA[1.0])

returns `f` and `∂f/∂T` at constant total volume and composition, where f is the total helmholtz energy, given by `eos(model,V,T,z)`

"""
function ∂f∂T(model,V,T,z=SA[1.0])
    return ForwardDiff.derivative(∂T -> eos(model,V,∂T,z),T)
end

"""
    ∂f∂V(model,V,T,z=SA[1.0])

returns `f` and `∂f/∂V` at constant temperature and composition, where f is the total helmholtz energy, given by `eos(model,V,T,z)`, and V is the total volume

"""
function ∂f∂V(model,V,T,z)
    return ForwardDiff.derivative(∂V -> eos(model,∂V,T,z),V)
end

#returns a tuple of the form ([∂f∂V,∂f∂T],f),using the least amount of computation
"""
    ∂f(model,V,T,z)

returns zeroth order (value) and first order derivative information of the total helmholtz energy (given by `eos(model,V,T,z)`).
the result is given in two values:

```
grad_f,fval = ∂2f(model,V,T,z)
```

where:
```
fval   = f(V,T) = eos(model,V,T,z)

grad_f = [ ∂f/∂V; ∂f/∂T]

 ```

Where `V` is the total volume, `T` is the temperature and `f` is the total helmholtz energy.
"""
function ∂f(model,V,T,z)
    f(w) = eos(model,first(w),last(w),z)
    V,T = promote(V,T)
    VT_vec = SVector(V,T)
    ∂result = DiffResults.GradientResult(VT_vec)
    res_∂f =  ForwardDiff.gradient!(∂result, f,VT_vec)
    _f =  DiffResults.value(res_∂f)
    _∂f = DiffResults.gradient(res_∂f)
    return (_∂f,_f)
end

#returns p and ∂p∂V at constant T
#it doesnt do a pass over temperature, so its
#faster that d2f when only requiring d2fdV2

"""
    p∂p∂V(model,V,T,z=SA[1.0])

returns `p` and `∂p/∂V` at constant temperature, where p is the pressure = `pressure(model,V,T,z)` and `V` is the total Volume.

"""
function p∂p∂V(model,V,T,z=SA[1.0])
    V_vec =   SVector(V)
    f(∂V) = pressure(model,only(∂V),T,z)
    ∂result = DiffResults.GradientResult(V_vec)
    res_∂f =  ForwardDiff.gradient!(∂result, f,V_vec)
    _p =  DiffResults.value(res_∂f)
    _∂p∂V = only(DiffResults.gradient(res_∂f))
    return _p,_∂p∂V
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
    f(w) = eos(model,first(w),last(w),z)
    V,T = promote(V,T)
    VT_vec =   SVector(V,T)
    ∂result = DiffResults.HessianResult(VT_vec)
    res_∂f =  ForwardDiff.hessian!(∂result, f,VT_vec)
    _f =  DiffResults.value(res_∂f)
    _∂f = DiffResults.gradient(res_∂f)
    _∂2f = DiffResults.hessian(res_∂f)
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
    f(w) = pressure(model,first(w),last(w),z)
    V,T = promote(V,T)
    VT_vec =   SVector(V,T)
    ∂result = DiffResults.HessianResult(VT_vec)
    res_∂f =  ForwardDiff.hessian!(∂result, f,VT_vec)
    _f =  DiffResults.value(res_∂f)
    _∂f = DiffResults.gradient(res_∂f)
    _∂2f = DiffResults.hessian(res_∂f)
    return (_∂2f,_∂f,_f)
end


"""
    f_hess(model,V,T,z)

returns the second order volume (`V`) and temperature (`T`) derivatives of the total helmholtz energy (given by `eos(model,V,T,z)`). the result is given in a 2x2 `SMatrix`, in the form:

```
[ ∂²f/∂V²  ∂²f/∂V∂T
 ∂²f/∂V∂T  ∂²f/∂V²]
 ```

use this instead of the ∂2f if you only need second order information. ∂2f also gives zeroth and first order derivative information, but due to a bug in the used AD, it allocates more than necessary.
"""
function f_hess(model,V,T,z)
    f(w) = eos(model,first(w),last(w),z)
    V,T = promote(V,T)
    VT_vec = SVector(V,T)
    return ForwardDiff.hessian(f,VT_vec)
end

"""
    ∂²³f(model,V,T,z=SA[1.0])

returns `∂²A/∂V²` and `∂³A/∂V³`, in a single ForwardDiff pass. used mainly in `crit_pure` objective function

"""
function ∂²³f(model,V,T,z=SA[1.0])
    V_vec =   SVector(V)
    f(∂A∂V) = ForwardDiff.derivative(∂²A∂V² -> pressure(model,only(∂²A∂V²),T,z),only(∂A∂V))
    ∂result = DiffResults.GradientResult(V_vec)
    res_∂f =  ForwardDiff.gradient!(∂result,f,V_vec)
    _∂²A∂V² =  DiffResults.value(res_∂f)
    _∂³A∂V³ = only(DiffResults.gradient(res_∂f))
    return _∂²A∂V², _∂³A∂V³
end
