"""
    SRK(components::Vector{String};
    idealmodel = BasicIdeal,
    alpha = SoaveAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[], 
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

## Description
Soave-Redlich-Kwong equation of state. it uses the following models:

- Translation Model: [`NoTranslation`](@ref)
- Alpha Model: [`SoaveAlpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

## Model Construction Examples
```julia
# Using the default database
model = SRK("water") #single input
model = SRK(["water","ethanol"]) #multiple components
model = SRK(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model
model = SRK(["water","ethanol"],alpha = Soave2019) #modifying alpha function, using 2019 correlation
model = SRK(["water","ethanol"],translation = RackettTranslation) #modifying translation
model = SRK(["water","ethanol"],mixing = KayRule) #using another mixing rule
model = SRK(["water","ethanol"],mixing = WSRule, activity = NRTL) #using advanced EoS+gᴱ mixing rule

# Passing a prebuilt model

my_alpha = SoaveAlpha(["ethane","butane"],userlocations = Dict(:acentricfactor => [0.1,0.2]))
model =  SRK(["ethane","butane"],alpha = my_alpha)

# User-provided parameters, passing files or folders

model = SRK(["neon","hydrogen"]; userlocations = ["path/to/my/db","cubic/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = SRK(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21]
                        k = [0. 0.18; 0.18 0.], #k,l can be ommited in single-component models.
                        l = [0. 0.01; 0.01 0.])
                    )
```


## References
1. Soave, G. (1972). Equilibrium constants from a modified Redlich-Kwong equation of state. Chemical Engineering Science, 27(6), 1197–1203. [doi:10.1016/0009-2509(72)80096-4](https://doi.org/10.1016/0009-2509(72)80096-4)
"""
function SRK(components;
    idealmodel = BasicIdeal,
    alpha = SoaveAlpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[], 
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    return RK(components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing = mixing,
    activity = activity,
    translation = translation,
    userlocations = userlocations, 
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    activity_userlocations = activity_userlocations,
    translation_userlocations = translation_userlocations,
    reference_state = reference_state,
    verbose = verbose)
end
export SRK