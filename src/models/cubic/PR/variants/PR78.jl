"""
    PR78(components::Vector{String};
    idealmodel = BasicIdeal,
    alpha = PR78Alpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[], 
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

Peng Robinson (1978) equation of state. It uses the following models:
- Translation Model: [`NoTranslation`](@ref)
- Alpha Model: [`PR78Alpha`](@ref)
- Mixing Rule Model: [`vdW1fRule`](@ref)

## Model Construction Examples
```julia
# Using the default database
model = PR78("water") #single input
model = PR78(["water","ethanol"]) #multiple components
model = PR78(["water","ethanol"], idealmodel = ReidIdeal) #modifying ideal model
model = PR78(["water","ethanol"],alpha = Soave2019) #modifying alpha function
model = PR78(["water","ethanol"],translation = RackettTranslation) #modifying translation
model = PR78(["water","ethanol"],mixing = KayRule) #using another mixing rule
model = PR78(["water","ethanol"],mixing = WSRule, activity = NRTL) #using advanced EoS+gᴱ mixing rule

# Passing a prebuilt model

my_alpha = PRAlpha(["ethane","butane"],userlocations = Dict(:acentricfactor => [0.1,0.2]))
model = PR78(["ethane","butane"],alpha = my_alpha) #this model becomes a normal PR EoS

# User-provided parameters, passing files or folders
model = PR78(["neon","hydrogen"]; userlocations = ["path/to/my/db","cubic/my_k_values.csv"])

# User-provided parameters, passing parameters directly

model = PR78(["neon","hydrogen"];
        userlocations = (;Tc = [44.492,33.19],
                        Pc = [2679000, 1296400],
                        Mw = [20.17, 2.],
                        acentricfactor = [-0.03,-0.21]
                        k = [0. 0.18; 0.18 0.], #k,l can be ommited in single-component models.
                        l = [0. 0.01; 0.01 0.])
                    )
```


## References
1. Robinson DB, Peng DY. The characterization of the heptanes and heavier fractions for the GPA Peng-Robinson programs. Tulsa: Gas Processors Association; 1978
"""
function PR78(components;
    idealmodel = BasicIdeal,
    alpha = PR78Alpha,
    mixing = vdW1fRule,
    activity = nothing,
    translation = NoTranslation,
    userlocations = String[], 
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    return PR(components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing = mixing,
    activity = activity,
    translation = translation,
    userlocations = userlocations,
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    translation_userlocations = translation_userlocations,
    reference_state = reference_state,
    verbose = verbose)
end
export PR78