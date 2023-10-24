# v0.5.6

## New Features
- New EoS: Polar PCSAFT with Quadrupolar interactions (`QPPCSAFT`)
- `Unitful.jl` support has been moved into an extension.
- `MultiComponentFlash.jl` support via extension: you can pass `model::EoSModel` to `MultiComponentFlash.flash_2ph`. there is also a new `MCFlashJL` method that calls `flash_2ph` using the `Clapeyron.tp_flash` interface.
- custom types can be passed to the `userlocations` keyword argument, defining `Clapeyron.can_nt(::datatype) = true` and `Clapeyron.to_nt(x::datatype)::Union{AbstractDict,NamedTuple}`

## Bug Fixes
- bug in `pharmaPCSAFT` mixing rules
