# v0.4.11

## New Features

- The package now performs precompilation (via `PrecompileTools.jl`) of some commonly used functions and EoSModels. this will add more time at the installation of the package, in exchange with decreased loading times each time `Clapeyron` is loaded in a new julia session. you can turn on/off the feature with `Clapeyron.precompile_clapeyron!(setting::Bool)` (recomended in the case of developing the library). due to how precompilation is done, it is only done from julia 1.9 onwards. 
- New EoS: NRTL with aspen paranetrization of τᵢⱼ : `aspenNRTL`
- `split_model` should be a little bit faster, and will perform correct splitting of group association models that use `get_group_idx`.

## Bug fixes

- Combining rule for epsilon in `SAFTgammaMie` updated from `√(ϵᵢ*ϵⱼ)*(σᵢ^3 * σⱼ^3)/σᵢⱼ^6` to `√(ϵᵢ*ϵⱼ*(σᵢ^3 * σⱼ^3))/σᵢⱼ^3`. the old behaviour can be obtained by passing the argument `epsilon_mixing::Symbol` to the model constructor. (#179)
