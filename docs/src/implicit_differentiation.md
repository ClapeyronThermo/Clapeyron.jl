# Implicit differentiation of root solvers

## Implicit function theorem
Given an implicit variable $x$ defined through the root $F(x(\theta),\theta)=0$, one can apply the implicit function theorem (IFT) to obtain the gradients of $x$ wrt $\theta$

```math
\frac{d F}{d \theta} = \frac{\partial F}{\partial x} \frac{d x}{d \theta} + \frac{\partial F}{\partial \theta} = 0 \implies \frac{d x}{d \theta} = - \left(\frac{\partial F}{\partial x}\right)^{-1} \frac{\partial F}{\partial \theta}
```

To avoid the unrolling of iterations, computing the derivative at each iteration and hence populating the AD stack unnecessarily, `Clapeyron` uses [`IFTDual.jl`](https://github.com/Garren-H/IFTDuals.jl), which implements the IFT for dual numbers and is hence only limited to `ForwardDiff.jl`. `IFTDuals.jl` can and is used to obtain abitrary order derivatives.

`Clapeyron` previously implemented a final Newton step to propagate first order derivatives, which is also limited to first order AD. Although this method does produce numeric values for nested Duals, these higher order derivatives are incorrect as it neglects the indirect effects from lower order gradients. For a more complete explanation see [this discussion](https://github.com/ClapeyronThermo/Clapeyron.jl/issues/502).
