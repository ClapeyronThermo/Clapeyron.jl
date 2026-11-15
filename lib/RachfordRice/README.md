# RachfordRice.jl

Rachford-Rice solver in julia.

## Quickstart

```julia
using RachfordRice
using RachfordRice: residual, rr_limits #not exported
k1 = [3.0E+00, 2.0E+00, 1.1E+00, 8.0E-01, 4.0E-01, 1.0E-01, 5.0E-02, 2.0E-02, 1.0E-02, 1.0E-04]
z1 = [0.30,0.15, 0.05, 0.02, 0.01, 0.02, 0.02, 0.03, 0.07, 0.33]

rr = RachfordRiceResidual(k1,z1) #generate a rachford rice residual object
k_i1 = rr[1] #indexing returns the inner K-value
Vlo,Vhi = rr_limits(rr) #Whitson-Michelsen bounds, if they exist.

g0 = residual(rr,0)
g1 = residual(rr,1)
result = rachfordrice(rr)
result2 = rachfordrice(rr,negative_flash = true) #support for negative flash
V = result.V #vapour fraction
L = result.L #liquid fraction
status = result.status #check if we have a liquid, a vapour, equilibria, a trivial case, or a failure.
feasible = result.feasible #if we have a single phase, check if there is any feasible solution to the rachford-rice problem, even outside (0,1), also known as negative flash.
x = liquid_compositions(rr,V)
x = liquid_compositions!(x,rr,V)
y = vapour_compositions(rr,V)
y = vapour_compositions!(y,rr,V)
```

## Features

- Support for "masking": rachford-rice is tipically used in the inner loop of a phase equilibria clculation, where the function is called multiple times with different values of `K` and `z`. on some cases, we know that one particular value of `K` is zero (a non-volatile component, it does not appear on `y`) or infinity (a non-condensable component, it does not appear on `x`). `RachfordRiceResidual` supports passing `non_inx` and `non_iny` keywords to allow specific components to be treated as if the K value at that particular index is zero or infinity.

  ```julia
  #mark a component as non-volatile
  non_volatile = fill(false,10)
  non_volatile[1] = true
  rr = RachfordRiceResidual(x,y;non_iny = non_volatile)
  rr[1] == 0 #true
  ``` 

- Exact solver for up to 4 components: The n-component Rachford-Rice problem can be reformulated as the root of a n-1 degree polynomial. `RachfordRice` takes advantage of this fact to solve the problem exactly when those conditions are met. Particularly,multiple non-volatile and non-condensable components are all lumped into one for each type.

## General Solver

If the rachford-rice problem is too big to be solved analitically, then the solution is obtained iteratively. We use a guarded halley approach: starting from the midpoint of the Whitson-Michelsen bounds (or an user supplied value), the next point is approximated via Halley; if the next iterate falls outside the bounds, then bisection is used as fallback.
