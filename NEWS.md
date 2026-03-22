# v0.6.23

## New Features

- New ideal models: `GCAlyLee` and `BurkhardtIdeal`
-New string macros: `@cas_str` and `@smiles_str`, that allow searching components by CAS and SMILES respectively.
- improved initial point for bubble/dew calculations

## Bug fixes

- Fixes and improvements in volume calculation.
- `SingleFluid`, when passing a JSON directly as a string, now the name stored is extracted and used. before, the whole JSON string was used as a name.
