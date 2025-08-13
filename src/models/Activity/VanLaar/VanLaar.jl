export VanLaar

struct VanLaarParam <: EoSParam
    A12::SingleParam{Float64}
    A21::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type VanLaarModel <: ActivityModel end

struct VanLaar{C<:EoSModel} <: VanLaarModel
    components::Vector{String}
    params::VanLaarParam
    puremodel::EoSVectorParam{C}
    references::Vector{String}
end
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
default_locations(::Type{VanLaar}) = [
    "properties/critical.csv",
    "properties/molarmass.csv",
    "Activity/VanLaar/vanlaar_unlike.csv"
]

function VanLaar(components;
    puremodel=PR,
    userlocations=String[],
    pure_userlocations=String[],
    verbose=false,
    reference_state=nothing)

    formatted_components = format_components(components)
    @assert length(formatted_components)==2

    params = getparams(
        formatted_components,
        default_locations(VanLaar);
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
    vp = VanLaarParam(A12,A21,Mw)
    v  = VanLaar(formatted_components, vp, pure, String[""])
    set_reference_state!(v,reference_state;verbose=verbose)
    binary_component_check(VanLaar,v)
    return v
end

function excess_g_vanlaar(m::VanLaarModel, p, T, z)
    n = sum(z)
    x = z ./ n
    A12 = m.params.A12.values[1]
    A21 = m.params.A21.values[1]
    ge = (A12*A21*x[1]*x[2]) / (A12*x[1] + A21*x[2])
    return n*R̄*T*ge
end

excess_gibbs_free_energy(model::VanLaarModel,p,T,z) = excess_g_vanlaar(model,p,T,z)
