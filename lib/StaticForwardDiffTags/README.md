# StaticForwardDiffTags

This package is composed of three related utilities

- `WithContext`: a wrapper over a callable struct that allows ForwardDiff to actually know whats inside
- `STag{}`: a variant of `ForwardDiff.Tag` that can be compared at compilation time, only generated for `WithContext` callables.
- `λFn` a simple function type wrapper to mark that functions are "pure", in the sense that they have no inaccessible inner state.

## Usage (proposed)

```julia
using StaticForwardDiffTags, ForwardDiff
f(p) = x -> sin(p*x) + p
fc = auto_context(f) #uses a generated function to recurse over the closure and obtain the number type of the captured state.

dfdx(x) = ForwardDiff.derivative(fc(2.0),x)
```
