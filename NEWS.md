# v0.4.11

## Bug fixes

- Combining rule for epsilon in `SAFTgammaMie` updated from `√(ϵᵢ*ϵⱼ)*(σᵢ^3 * σⱼ^3)/σᵢⱼ^6` to `√(ϵᵢ*ϵⱼ*(σᵢ^3 * σⱼ^3))/σᵢⱼ^3`. the old behaviour can be obtained by passing the argument `epsilon_mixing::Symbol` to the model constructor. (#179)
