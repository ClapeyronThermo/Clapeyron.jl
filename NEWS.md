# v0.6.22

## New Features

- Implicit differentiation: the core implicit differentiation routines were replaced by `IFTDuals.ift`.
- Experimental: `AssocOptions` supports the `implicit_ad` option to calculate derivatives of the association solver via `IFTDuals`
- Implicit differentiation is now enabled in PH flash and PS flash.
- Experimental: new tpd function: `Clapeyron.tpd2`, that returns a `TPDResult` struct instead of a tuple of vectors
- Experimental: New model wrapper for electrolyte wrappers: `MeanIonicApproach` with support for `tp_flash`

## Bug fixes

- Fixes in `iPCSAFT`
- Fixes in electrolyte routines
- Fixes in Multifluid initial volume
- Fixes in `MultiPhaseTPFlash`
