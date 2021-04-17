#derivative logic

#ForwardDiff compiles one method per function, so separating the
#differentiation logic from the property logic allows the differentials
#to be compiled only once


"""
    ∂f∂t(model,v,t,z=SA[1.0])

returns `f` and `∂f/∂T` at constant total volume and composition, where f is the total helmholtz energy, given by `eos(model,v,T,z)`

"""
function ∂f∂t(model,v,t,z=SA[1.0])
    return ForwardDiff.derivative(∂t -> eos(model,v,∂t,z),t)
end

"""
    ∂f∂v(model,v,t,z=SA[1.0])

returns `f` and `∂f/∂v` at constant temperature and composition, where f is the total helmholtz energy, given by `eos(model,v,T,z)`, and v is the total volume

"""
function ∂f∂v(model,v,t,z)
    return ForwardDiff.derivative(∂v -> eos(model,∂v,t,z),v)
end

#returns a tuple of the form ([∂f∂v,∂f∂t],f),using the least amount of computation
"""
    ∂f(model,v,t,z)

returns zeroth order (value) and first order derivative information of the total helmholtz energy (given by `eos(model,v,T,z)`).
the result is given in two values:

```
grad_f,fval = ∂2f(model,v,t,z)
```

where:
```
fval   = f(v,T) = eos(model,v,T,z)

grad_f = [ ∂f/∂v; ∂f/∂T]

 ```

Where `v` is the total volume, `T` is the temperature and `f` is the total helmholtz energy.
"""
function ∂f(model,v,t,z)
    f(w) = eos(model,first(w),last(w),z)
    v,t = promote(v,t)
    vt_vec = SVector(v,t)
    ∂result = DiffResults.GradientResult(vt_vec)
    res_∂f =  ForwardDiff.gradient!(∂result, f,vt_vec)
    _f =  DiffResults.value(res_∂f)
    _∂f = DiffResults.gradient(res_∂f)
    return (_∂f,_f)
end

#returns p and ∂p∂v at constant t
#it doesnt do a pass over temperature, so its
#faster that d2f when only requiring d2fdv2

"""
    p∂p∂v(model,v,t,z=SA[1.0])

returns `p` and `∂p/∂v` at constant temperature, where p is the pressure = `pressure(model,v,T,z)` and `v` is the total volume.

"""
function p∂p∂v(model,v,t,z=SA[1.0])
    v_vec =   SVector(v)
    f(∂v) = pressure(model,only(∂v),t,z)
    ∂result = DiffResults.GradientResult(v_vec)
    res_∂f =  ForwardDiff.gradient!(∂result, f,v_vec)
    _p =  DiffResults.value(res_∂f)
    _∂p∂v = only(DiffResults.gradient(res_∂f))
    return _p,_∂p∂v
end

"""
    ∂2f(model,v,t,z)

returns zeroth order (value), first order and second order derivative information of the total helmholtz energy (given by `eos(model,v,T,z)`).
the result is given in three values:

```
hess_f,grad_f,fval = ∂2f(model,v,t,z)
```

where:
```
fval   = f(v,T) = eos(model,v,T,z)

grad_f = [ ∂f/∂v; ∂f/∂T]

hess_f = [ ∂²f/∂v²; ∂²f/∂v∂T
          ∂²f/∂v∂T; ∂²f/∂v²]
 ```

Where `v` is the total volume, `T` is the temperature and `f` is the total helmholtz energy.
"""
function ∂2f(model,v,t,z)
    f(w) = eos(model,first(w),last(w),z)
    v,t = promote(v,t)
    vt_vec =   SVector(v,t)
    ∂result = DiffResults.HessianResult(vt_vec)
    res_∂f =  ForwardDiff.hessian!(∂result, f,vt_vec)
    _f =  DiffResults.value(res_∂f)
    _∂f = DiffResults.gradient(res_∂f)
    _∂2f = DiffResults.hessian(res_∂f)
    return (_∂2f,_∂f,_f)
end

"""
    ∂2p(model,v,t,z)

returns zeroth order (value), first order and second order derivative information of the pressure.
the result is given in three values:

```
hess_p,grad_p,pval = ∂2p(model,v,t,z)
```

where:
```
pval   = p(v,T) = pressure(model,v,T,z)

grad_p = [ ∂p/∂v; ∂p/∂T]

hess_p = [ ∂²p/∂v²; ∂²p/∂v∂T
          ∂²p/∂v∂T; ∂²p/∂v²]
 ```

Where `v` is the total volume, `T` is the temperature and `p` is the pressure.
"""
function ∂2p(model,v,t,z)
    f(w) = pressure(model,first(w),last(w),z)
    v,t = promote(v,t)
    vt_vec =   SVector(v,t)
    ∂result = DiffResults.HessianResult(vt_vec)
    res_∂f =  ForwardDiff.hessian!(∂result, f,vt_vec)
    _f =  DiffResults.value(res_∂f)
    _∂f = DiffResults.gradient(res_∂f)
    _∂2f = DiffResults.hessian(res_∂f)
    return (_∂2f,_∂f,_f)
end


"""
    f_hess(model,v,t,z)

returns the second order volume (`v`) and temperature (`T`) derivatives of the total helmholtz energy (given by `eos(model,v,T,z)`). the result is given in a 2x2 `SMatrix`, in the form:

```
[ ∂²f/∂v²  ∂²f/∂v∂T
 ∂²f/∂v∂T  ∂²f/∂v²]
 ```

use this instead of the ∂2f if you only need second order information. ∂2f also gives zeroth and first order derivative information, but due to a bug in the used AD, it allocates more than necessary.
"""
function f_hess(model,v,t,z)
    f(w) = eos(model,first(w),last(w),z)
    v,t = promote(v,t)
    vt_vec = SVector(v,t)
    return ForwardDiff.hessian(f,vt_vec)
end

"""
    ∂p2∂p3(model,v,t,z=SA[1.0])

returns `∂²p/∂v²` and `∂³p/∂v³`, in a single ForwardDiff pass. used mainly in `crit_pure` objective function

"""
function ∂p2∂p3(model,v,t,z=SA[1.0])
    v_vec =   SVector(v)
    f(∂v) = ForwardDiff.derivative(∂2v -> pressure(model,only(∂2v),t,z),only(∂v))
    ∂result = DiffResults.GradientResult(v_vec)
    res_∂f =  ForwardDiff.gradient!(∂result, f,v_vec)
    _p =  DiffResults.value(res_∂f)
    _∂p∂v = only(DiffResults.gradient(res_∂f))
    return _p,_∂p∂v
end
