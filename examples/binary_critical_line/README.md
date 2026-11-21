# Trace a local binary critical line

Run from the Clapeyron.jl checkout with its project dependencies instantiated:

```sh
julia --project=. --compile=min -O0 examples/binary_critical_line/run.jl
julia --project=. --compile=min -O0 examples/binary_critical_line/check.jl
```

The example starts from a known interior PC-SAFT methane/decane root and takes
short pseudo-arclength steps in both directions using log temperature, log
density and mole fraction. A continuation curve can turn in composition without
forcing the solver to use composition as its independent variable.

The corrector uses a composition-scaled Helmholtz Hessian. When the transverse
eigenvalue is positive and separated from zero, dividing the two residuals by
that eigenvalue and its square leaves their zero set unchanged. The result
reports raw and normalized values separately; using the same numerical
tolerance does not imply the same stopping behavior.

The script returns two local critical endpoints near 578.71 and 579.87 K on
its selected temperature boundaries. This is an example with an explicit seed
and model inputs, not a global phase-stability test, a promise to find all
branches, or a supported API. Only the selected coordinate is corrected to an
endpoint. Degenerate critical modes and pure-component endpoints are outside
the example's scope. The checked values are model calculations, not experimental
data or a fit-range validation. No timing comparison is implied by the compiler
flags.

The fixed-composition equations are adapted from Clapeyron's MIT-licensed
implementation. See the [project license](../../LICENSE.md). For the wider
critical-line question and literature references, see [#318](https://github.com/ClapeyronThermo/Clapeyron.jl/issues/318).
