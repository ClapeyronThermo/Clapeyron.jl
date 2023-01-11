# v0.4.0

## New Features

- New estimation framework, to fit new model parameters from known data. at the moment of this release, it should support all SAFTs, Cubics, Activity models and models that don't require any special pre-computed cache.

- `Base.transpose(model::EoSModel)` is defined. now this is valid code:
```julia
T = 300:350
comp = ["water"]
models = [m(comp) for m in (PR,PCSAFT,SAFTVRMie)]
saturation_pressure.(transpose(models),T) #produces a 51×3 Matrix{Tuple{Float64, Float64, Float64}}:
```

## Breaking Changes
- `x0_sat_pure` now returns `Vl` and `Vv` (in `[m^3]`) instead of `log10(Vl)` and `log10(Vv)`.
- in databases, `segment` is used instead of `m`
- in databases, `Vc` is used instead of `vc`
- in databases, `Pc` is used instead of `pc`
- in databases, `acentricfactor` is used instead of `w`


