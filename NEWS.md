# v0.6.24

## New Features

- Association: better initial points and faster evaluation via compression of the association matrix.
- `GeneralizedXYFlash`: added `verbose` option
- `ChemPotBubblePressure`: added `verbose` option
- `ChemPotBubbleTemperature`: added `verbose` option
- `ChemPotDewPressure`: added `verbose` option
- `ChemPotDewTemperature`: added `verbose` option
- `Tproperty`/`Pproperty`: better initial points for poins inside the phase change region.
- `spinodal_pressure`/`spinodal_temperature`: initial points now use `edge_temperature`/`edge_pressure` instead of bubble/dew calculations, improving speed and stability, especially with conditions near the mixture critical point
- `EPPR78`: support for passing `alpha` and `translation` submodels
- stability improvements to `edge_pressure`/`edge_temperature`

## Bug fixes

- Fixes in `split_model` when indices aren't ordered
