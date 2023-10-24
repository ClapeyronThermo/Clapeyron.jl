# v0.5.7

## New Features
- New EoS: Polar PCSAFT with Quadrupolar interactions (`QPPCSAFT`)
- New EoS: Group-Contribution simplified PC-SAFT (`gcsPCSAFT`)
- New EoS: Group-Contribution homosegmented polar PC-SAFT (`gcPPCSAFT`)
- `Unitful.jl` support has been moved into an extension.
- `MultiComponentFlash.jl` support via extension: you can pass `model::EoSModel` to `MultiComponentFlash.flash_2ph`. there is also a new `MCFlashJL` method that calls `flash_2ph` using the `Clapeyron.tp_flash` interface.
- custom types can be passed to the `userlocations` keyword argument, defining `Clapeyron.can_nt(::datatype) = true` and `Clapeyron.to_nt(x::datatype)::Union{AbstractDict,NamedTuple}`

## Bug Fixes
- bug in `pharmaPCSAFT` mixing rules
