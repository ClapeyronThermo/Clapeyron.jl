struct VanLaarParam <: EoSParam
    A₁₂::SingleParam{Float64}
    A₂₁::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type VanLaarModel <: ActivityModel end

struct VanLaar{c<:EoSModel} <: VanLaarModel
    components::Array{String,1}
    params::VanLaarParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export VanLaar

"""
    VanLaar <: ActivityModel
    VanLaar(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)

## Input parameters
- `A12`: Single Parameter (`Float64`, defaults to `0`) - Binary Interaction Parameter
- `A21`: Single Parameter (`Float64`, defaults to `0`) - Binary Interaction Parameter
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`

## model parameters
- `A12`: Single Parameter (`Float64`, defaults to `0`) - Binary Interaction Parameter
- `A21`: Single Parameter (`Float64`, defaults to `0`) - Binary Interaction Parameter

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
VanLaar activity coefficient model, for binary mixture:
```
Gᴱ = nRT·(A12·x₁·A21·x₂/(A12·x₁ + A21·x₂))
```

## Model Construction Examples
```
# Using the default database
model = VanLaar(["water","ethanol"]) #Default pure model: PR
model = VanLaar(["water","ethanol"],puremodel = BasicIdeal) #Using Ideal Gas for pure model properties
model = VanLaar(["water","ethanol"],puremodel = PCSAFT) #Using Real Gas model for pure model properties

# Passing a prebuilt model

my_puremodel = AbbottVirial(["water","ethanol"]; userlocations = ["path/to/my/db","critical.csv"])
mixing = VanLaar(["water","ethanol"],puremodel = my_puremodel)

# Using user-provided parameters

# Passing files or folders
model = VanLaar(["water","ethanol"];userlocations = ["path/to/my/db","vanLaar.csv"])

# Passing parameters directly
model = VanLaar(["water","ethanol"],
        userlocations = (A12 = [4512],
                        A21 = [3988.52],
                        Mw = [18.015, 46.069])
                        )
```

## References
[1] J. J. van Laar, Sechs Vorträgen über das thermodynamische Potential. (Six Lectures on the Thermodynamic Potential). Braunschweig, Fried. Vieweg & Sohn, 1906.
[2] J. J. van Laar, “Über Dampfspannungen von binären Gemischen (The vapor pressure of binary mixtures)”, Z. Physik. Chem., vol. 72, pp. 723-751, May 1910.
[3] J. J. van Laar, “Zur Theorie der Dampfspannungen von binären Gemischen. Erwiderung an Herrn F. Dolezalek (Theory of vapor pressure of binary mixtures. Reply to Mr. F. Dolezalek)”, Z. Physik. Chem., vol. 83, pp. 599-608, June 1913. 
"""
VanLaar

default_locations(::Type{VanLaar}) = ["properties/critical.csv", "properties/molarmass.csv","Activity/VanLaar/vanLaar_unlike.csv"]

function VanLaar(components;
    puremodel = PR,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false,
    reference_state = nothing)
    
    formatted_components = format_components(components)
    params = getparams(formatted_components, default_locations(VanLaar); userlocations = userlocations, asymmetricparams=["A12","A21"], ignore_missing_singleparams=["A12","A21","Mw"], verbose = verbose)
    A12        = params["A12"]
    A21        = params["A21"]
    Mw         = params["Mw"]
    
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    packagedparams = VanLaarParam(A12,A21,Mw)
    references = String[""]
    model = VanLaar(formatted_components,packagedparams,_puremodel,references)
    set_reference_state!(model,reference_state,verbose = verbose)
    binary_component_check(VanLaar,model)
    return model
end

function excess_g_VanLaar(model::VanLaarModel, p, T, z)
    n = sum(z)
    x = z ./ n

    A12 = model.params.A12.values[1]
    A21 = model.params.A21.values[1]

    ge = A12 * x[1] * A21 * x[2] / (A12*x[1] + A21*x[2])

    return n * R̄ * T * ge
end