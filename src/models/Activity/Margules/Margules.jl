export Margules

struct MargulesParam <: EoSParam
    A12::SingleParam{Float64}
    A21::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type MargulesModel <: ActivityModel end

struct Margules{C<:EoSModel} <: MargulesModel
    components::Vector{String}
    params::MargulesParam
    puremodel::EoSVectorParam{C}
    references::Vector{String}
end
"""
    Margules <: ActivityModel
    Margules(components;
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
Margules activity coefficient model, for binary mixture:
```
Gᴱ = nRT·(x₁·x₂·(A21·x₁+A12·x₂))
```
## Model Construction Examples
```
# Using the default database
model = Margules(["water","ethanol"]) #Default pure model: PR
model = Margules(["water","ethanol"],puremodel = BasicIdeal) #Using Ideal Gas for pure model properties
model = Margules(["water","ethanol"],puremodel = PCSAFT) #Using Real Gas model for pure model properties
# Passing a prebuilt model
my_puremodel = AbbottVirial(["water","ethanol"]; userlocations = ["path/to/my/db","critical.csv"])
mixing = Margules(["water","ethanol"],puremodel = my_puremodel)
# Using user-provided parameters
# Passing files or folders
model = Margules(["water","ethanol"];userlocations = ["path/to/my/db","margules.csv"])
# Passing parameters directly
model = Margules(["water","ethanol"],
        userlocations = (A12 = [4512],
                        A21 = [3988.52],
                        Mw = [18.015, 46.069])
                        )
```
## References
1. Max Margules, « Über die Zusammensetzung der gesättigten Dämpfe von Misschungen », Sitzungsberichte der Kaiserliche Akadamie der Wissenschaften Wien Mathematisch-Naturwissenschaftliche Klasse II, vol. 104, 1895, p. 1243–1278
"""
Margules
default_locations(::Type{Margules}) = [
    "properties/critical.csv",
    "properties/molarmass.csv",
    "Activity/Margules/margules_unlike.csv"
]

function Margules(components;
    puremodel=PR,
    userlocations=String[],
    pure_userlocations=String[],
    verbose=false,
    reference_state=nothing)

    formatted_components = format_components(components)
    @assert length(formatted_components)==2

    params = getparams(
        formatted_components,
        default_locations(Margules);
        userlocations=userlocations,
        asymmetricparams=String[],
        ignore_missing_singleparams=["Mw"],
        verbose=verbose
    )

    _fill(name,comps,v) = (sp=SingleParam(name,comps); sp.values .= float(v); sp)

    function _as_single(name,comps,x; ij=(1,2))
        if x===nothing
            return SingleParam(name,comps)
        elseif x isa SingleParam
            return x
        elseif x isa PairParam
            return _fill(name,comps,x.values[ij...])
        else
            throw(ArgumentError("Cannot convert $(typeof(x)) to SingleParam"))
        end
    end

    Araw   = get(params,"A",   nothing)
    A12raw = get(params,"A12", nothing)
    A21raw = get(params,"A21", nothing)
    Mw     = get(params,"Mw",  SingleParam("Mw",formatted_components))

    A12 = _as_single("A12",formatted_components, A12raw===nothing ? Araw : A12raw; ij=(1,2))
    A21 = _as_single("A21",formatted_components, A21raw===nothing ? Araw : A21raw; ij=(2,1))

    pure = init_puremodel(puremodel,components,pure_userlocations,verbose)
    mp = MargulesParam(A12,A21,Mw)
    m  = Margules(formatted_components, mp, pure, String[""])
    set_reference_state!(m,reference_state;verbose=verbose)
    binary_component_check(Margules,m)
    return m
end

function excess_g_margules(m::MargulesModel, p, T, z)
    n = sum(z)
    x = z ./ n
    A12 = m.params.A12.values[1]
    A21 = m.params.A21.values[1]
    ge = x[1]*x[2]*(A21*x[1] + A12*x[2])
    return n*R̄*T*ge
end
