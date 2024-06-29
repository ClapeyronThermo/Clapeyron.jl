# v0.6.3

- Association solver is now faster for small association matrices.
- Michelsen TP-Flash: in case of valid K values but single phase rachford-rice, the procedure will assume bubble or dew point as a first iteration.
- `split_model` now works for `ClapeyronParam`,`Symbol`,`Number`,`AbstractString`,`Tuple`,`Missing` and `Nothing`. before those could only be splitted if inside an `EoSModel`.
- New function: `split_model_binaries`, that returns a list of all binary combinations of an n-component model.

## Bug fixes
- `SAFTgammaMie` fixes.