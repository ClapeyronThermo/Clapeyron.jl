# StaticForwardDiffTags

## Usage

```julia

using StaticForwardDiffTags, ForwardDiff
f(x,p) = sin(p*x) + p

#fixed_f(x) = f(p)(x)
fixed_f = SDiffFunction(x -> Base.Fix2(f,x),2.0)
ff(p) = SDiffFunction(x -> Base.Fix2(f,x),p)(3.0)

dfdx(x) = ForwardDiff.derivative(fixed_f,x)
dfdp(x) = ForwardDiff.derivative(ff,x)
d2f(x) = ForwardDiff.derivative(dfdx,x)
f2fdxdp(x) = ForwardDiff.derivative(dfdp,x)

fixed_f(1.0) #equivalent to sin(2*x)
dfdx(1.0) 
dfdp(1.0)
d2f(1.0)
f2fdxdp(1.0)

```
