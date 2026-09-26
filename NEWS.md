# v0.6.29

## New Features

- XY flashes: support for ideal models
- New model: Tillner-Roth-Friend for water-ammonia mixtures (`TillnerRothFriend`)

## Bug fixes

- `volume`: faster evaluation of gibbs energy difference
- Multiparameter: improved lower bound volume suggestion
- misc equilibria fixes
- `cPR`: fixed constructor function
- activity models: fixed excess functions
- `UNIFACFV`/`UNIFACFVPoly`: fixed inconsistency in the combinatorial term. Results may change between the old and new version, with the new version being the correct one.
