# Getting Started - Model Construction

All functions in Clapeyron revolve around an object we refer to as the `model`.
These models are intended to hold all of the information required to model a given system within a specific equation of state.
The absolute simplest of these is the `BasicIdeal` model, which models species as an ideal gas with only translational modes of motion.
It can be constructed simply using:

```julia
julia> model = BasicIdeal(["water"])
BasicIdeal()
```

This `model` is unique as it is the only model object that does not hold any information as the ideal gas is a universal model that does not depend on the chemical identity of the species.
Nevertheless, one can use this model to obtain properties such as the volume and isobaric heat capacity at a given temperature and pressure:

```julia
julia> volume(model,1e5,298.15)
0.02478957029602388

julia> isobaric_heat_capacity(model,1e5,298.15)
20.786156545383097
```

At this stage, we point out that all values within Clapeyron are in SI units, unless otherwise specified.

From here, we will now consider how `model`s are constructed in different equations of state.

!!! tip "Defining your own parameters"
    For guidance on how to define your own parameters, please examine the [User-Defined Parameters](./user_defined_parameters.md) section of the docs.

## Generic Models

Although `BasicIdeal` does provide the universal definition for the ideal gas model, it does fall short in one aspect: accounting for modes of motion beyond translation, such as rotations and vibrations.
Accounting for these contributions does involve providing chemical-specific parameters.
Thankfully, Clapeyron has a built-in databank of parameters for all supported equations of state.
As an example, we consider carbon-dioxide modelled by the `ReidIdeal` model, as opposed to the basic ideal model:

```julia
julia> model1 = BasicIdeal(["carbon dioxide"])
BasicIdeal()

julia> model2 = ReidIdeal(["carbon dioxide"])
ReidIdeal with 1 component:
 "carbon dioxide"
Contains parameters: coeffs
```

The difference between these models is best observed when one considers the isobaric heat capacity:

```julia
julia> isobaric_heat_capacity(model1,1e5,298.15)
20.786156545383097

julia> isobaric_heat_capacity(model2,1e5,298.15)
37.16203981139605
```

The difference in values can be attributed to the missing rotational and vibration contributions.
This difference is important to bear in mind when considering similar properties in other equations of state.

Nevertheless, due to the need for species-specific parameters, this `model2` contains much more information than for the `BasicIdeal` model:

* Component names:

  ```julia
  julia> model2.components
  1-element Vector{String}:
  "carbon dioxide"
  ```

* References for the equation of state:

  ```julia
  julia> model2.references
  "10.1063/1.3060771"
  ```

* The species-specific parameters:

  ```julia
  julia> model2.params
  Clapeyron.ReidIdealParam for ["carbon dioxide"] with 1 param:
  coeffs::SingleParam{NTuple{4, Float64}}

  julia> model2.params.coeffs
  SingleParam{NTuple{4, Float64}}("Reid Coefficients") with 1 component:
  "carbon dioxide" => (19.8, 0.0734, -5.6e-5, 1.72e-8)
  ```

One can also create a model for mixtures in a similar fashion:

```julia
julia> model3 = ReidIdeal(["carbon dioxide","methane"])
ReidIdeal with 2 components:
 "carbon dioxide"
 "methane"
Contains parameters: coeffs

julia> model3.params.coeffs
SingleParam{NTuple{4, Float64}}("Reid Coefficients") with 2 components:
 "carbon dioxide" => (19.8, 0.0734, -5.6e-5, 1.72e-8)
 "methane" => (19.25, 0.0521, 1.2e-5, -1.13e-8)
```

Ideal gas models are by far the simplest case to consider when building models.
Other equations of state have additional features which will be discussed next.
Of interest to general users may be the [group-contribution models](#group-contribution-models).

## Cubic Models

!!! tip "List of Models"
    A full list of cubic equations of state is available (see [Cubics](../eos/cubic.md)).

At the surface, cubic models are quite simple as well.
As an example, consider a mixture of methanol and benzene in Peng-Robinson (`PR`):

```julia
julia> model = PR(["methanol","benzene"])
PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule} with 2 components:
 "methanol"
 "benzene"
Contains parameters: a, b, Tc, Pc, Mw

julia> model.components
2-element Vector{String}:
 "methanol"
 "benzene"

julia> model.references
1-element Vector{String}:
 "10.1021/I160057A011"

julia> model.params.a
2×2 PairParam{Float64}(["methanol", "benzene"]) with values:
 1.02049  1.4453
 1.4453   2.04695
```

Note that, as we are dealing with a mixture, we need to include binary parameters (parameters that depend on two species).
This is stored as shown in the case of `a`.

Cubics are unique due to their modular nature (we can swap out pieces of the equation to hopefully make a more-accurate model).
This includes:

* **Alpha functions**, which improve the prediction the pure-component saturation curves:

  ```julia
  julia> model = PR(["methanol","benzene"];alpha=TwuAlpha)
  PR{BasicIdeal, TwuAlpha, NoTranslation, vdW1fRule} with 2 components:
  "methanol"
  "benzene"
  Contains parameters: a, b, Tc, Pc, Mw

  julia> model.alpha
  TwuAlpha with 2 components:
  "methanol"
  "benzene"
  Contains parameters: M, N, L
  ```

* **Volume translation methods**, which improve the prediction of the liquid volume:

  ```julia
  julia> model = PR(["methanol","benzene"];translation=PenelouxTranslation)
  PR{BasicIdeal, PRAlpha, PenelouxTranslation, vdW1fRule} with 2 components:
  "methanol"
  "benzene"
  Contains parameters: a, b, Tc, Pc, Mw

  julia> model.translation
  PenelouxTranslation with 2 components:
  "methanol"
  "benzene"
  Contains parameters: Vc, v_shift
  ```

* **Mixing rules**, which can improve the predicted phase behaviour.
  These come in two flavours:

  * Generic one-fluid mixing rules such as the van der Waals one-fluid mixing rule (the default in all cubics):

    ```julia
    julia> model = PR(["methanol","benzene"])
    PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule} with 2 components:
    "methanol"
    "benzene"
    Contains parameters: a, b, Tc, Pc, Mw

    julia> model.mixing
    vdW1fRule()
    ```

  * ``G_E``-mixing rules where an activity coefficient model is used in conjunction with the cubic equation of state:

    ```julia
    julia> model = PR(["methanol","benzene"];mixing=HVRule,activity=UNIFAC)
    PR{BasicIdeal, PRAlpha, NoTranslation, HVRule{UNIFAC{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}}}} with 2 components:
    "methanol"
    "benzene"
    Contains parameters: a, b, Tc, Pc, Mw

    julia> model.mixing
    HVRule{UNIFAC{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}}} with 2 components:
    "methanol"
    "benzene"

    julia> model.mixing.activity
    UNIFAC{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}} with 2 components:
    "methanol": "CH3OH" => 1
    "benzene": "ACH" => 6
    Group Type: UNIFACDortmund
    Contains parameters: A, B, C, R, Q
    ```

Whilst one could combine all the parts listed above in endless ways, there are some default combinations which we provide.
These are typically referred to as predictive cubics: Predictive SRK (`PSRK`) and Volume-Translated PR (`VTPR`).

## Activity Coefficient Models

!!! tip "List of Models"
    A full list of activity coefficient models is available (see [Activity Models](../eos/activity.md)).

As hinted at above, Clapeyron also supports activity coefficient models.
These can be assembled in a similar fashion to our previous models:

```julia
julia> model = Wilson(["water","ethanol"])
Wilson{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}} with 2 components:
 "water"
 "ethanol"
Contains parameters: g, Tc, Pc, ZRA, Mw

julia> model.params.g
2×2 PairParam{Float64}(["water", "ethanol"]) with values:
    0.0   3988.52
 1360.12     0.0

julia> model.puremodel
Clapeyron.EoSVectorParam{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}} with 2 components:
 "water"
 "ethanol"
```

The first noteworthy change is that, where most binary parameters we considered were previously symmetric, activity coefficient models have asymmetric parameters ($g_{ij}\neq g_{ji}$).
More importantly, while activity coefficient models are intended to model the liquid phase, through Raoult's law, they can be used to model vapour–liquid equilibrium.
For this to work, they need a model which can predict the saturation pressure of the pure components.
This is why you'll field the `puremodel` field within the `model` for activity coefficient models (as shown above).
By default, this is set to Peng–Robinson.
However, one can change this using the `puremodel` optional argument:

```julia
julia> model = Wilson(["water","ethanol"]; puremodel=PCSAFT)
Wilson{PCSAFT{BasicIdeal, Float64}} with 2 components:
 "water"
 "ethanol"
Contains parameters: g, Tc, Pc, ZRA, Mw
```

This pure model plays an important role in modelling the bulk properties of activity coefficient models as these approaches only predict excess properties:
$$X^E = X^\text{mix.}-\sum_ix_i X^\text{pure}$$
To obtain $X^\text{mix.}$, $X^\text{pure}$ is obtained using the `puremodel`.
Note that this means that, for all activity coefficient models, as they are pressure/volume independent, all assume ideal volume of mixing.

!!! warning "Coming changes"
    In a future update, Activity Coefficient Models will only be able to model the liquid phase by default.
    To model both the vapour and liquid phase, users will need to construct the model using [Composite Models](../eos/correlations.md).

### COSMO-SAC Models

Clapeyron.jl also supports COSMO-SAC-based models.
However, we only provide the activity coefficient model and not the quantum chemistry-level calculations required to obtain the sigma profiles.
As such, the required parameters for COSMO-SAC are not the sigma profiles, which are stored as vectors:

```julia
julia> model = COSMOSAC02(["water","ethanol"])
COSMOSAC02{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}} with 2 components:
 "water"
 "ethanol"
Contains parameters: Pi, V, A

julia> model.params.Pi
SingleParam{Vector{Float64}}("Pi") with 2 components:
 "water" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.104369  …  2.153543, 0.524173, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 "ethanol" => [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.294716  …  0.755815, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```

Note that, in comparison to other activity coefficient models, COSMO-SAC models will be quite a bit slower.
Furthermore, due to the size of the sigma profiles, we do not store a full database of parameters locally.
The parameters are usually obtained from the [NIST database](https://github.com/usnistgov/COSMOSAC) by specifying the `use_nist_database=true` optional argument.
Please verify the NIST database's license before usage.

## SAFT Models

!!! tip "List of Models"
    A full list of SAFT equations of state is available (see [SAFT and CPA Models](../eos/saft.md)).

Where most models use either single or binary parameters, SAFT models introduce a third type of parameter: association parameters.
These can be best seen when building the model:

```julia
julia> model = PCSAFT(["water","1-propanol"])
PCSAFT{BasicIdeal, Float64} with 2 components:
 "water"
 "1-propanol"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> model.sites
SiteParam with 2 components:
 "water": "e" => 1, "H" => 1
 "1-propanol": "e" => 1, "H" => 1

julia> model.params.epsilon_assoc
AssocParam{Float64}["water", "1-propanol"]) with 2 values:
("water", "e") >=< ("water", "H"): 2500.7
("1-propanol", "e") >=< ("1-propanol", "H"): 2276.8
```

When adding association, we allow for interactions between sites on species, hence why the `sites` field has been added.
The association parameters are also stored in a slightly different way as well.

One thing to note is that, unless cross-associating parameters are provided within the database, Clapeyron will not assume any other association interactions.
To include these, one can use the `AssocOptions` struct where, as well as other numerical settings, users can specify to use a combining rule:

```julia
julia> model = PCSAFT(["water","1-propanol"];assoc_options=AssocOptions(combining=:esd))
PCSAFT{BasicIdeal, Float64} with 2 components:
 "water"
 "1-propanol"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> model.params.epsilon_assoc
AssocParam{Float64}["water", "1-propanol"]) with 4 values:
("water", "e") >=< ("water", "H"): 2500.7
("1-propanol", "e") >=< ("water", "H"): 2388.75
("1-propanol", "H") >=< ("water", "e"): 2388.75
("1-propanol", "e") >=< ("1-propanol", "H"): 2276.8
```

This may be important when trying to obtain accurate predictions for mixtures.

### Cubic Plus Association

Another class of equation of state that doesn't really fit in either cubics or SAFT equations is CPA, where the two approaches are combined.
As a result, the features available in both approaches are extended to CPA:

```julia
julia> model = CPA(["methanol","ethane"])
CPA{BasicIdeal, RK{BasicIdeal, CPAAlpha, NoTranslation, vdW1fRule}} with 2 components:
 "methanol"
 "ethane"
Contains parameters: a, b, c1, Tc, epsilon_assoc, bondvol, Mw

julia> model = CPA(["methanol","ethane"];cubicmodel=PR,mixing=HVRule,activity=UNIFAC)
CPA{BasicIdeal, PR{BasicIdeal, CPAAlpha, NoTranslation, HVRule{UNIFAC{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}}}}} with 2 components:
 "methanol"
 "ethane"
Contains parameters: a, b, c1, Tc, epsilon_assoc, bondvol, Mw
```

Making our CPA implementation one of the most-extensible available.

## Empirical Equations of State

Clapeyron also supports high-accuracy, empirical equations of state.
These differ from the previous equations of state primarily because of large number of parameters needed to model a single species.
An example of this would be IAPWS-95 (for water):

```julia
julia> model = IAPWS95()
MultiParameter Equation of state for water:
 Polynomial power terms: 7
 Exponential terms: 44
 Gaussian bell-shaped terms: 3
 Non Analytic terms: 2
```

As these equations of state typically have common terms, rather than specifying the parameters, we highlight what terms the equation of state is made up of, as shown above.

!!! tip "List of Models"
    A full list of Empirical equations of state is available (see [Empirical Helmholtz Models](../eos/empiric.md)).
    The list of available systems can be expanded by including the CoolProp.jl extension (see [Extension - CoolProp](./extensions_coolprop.md)).

## Composite Models

Not every equation of state provides a global representation of the phase space of a system (for example, activity coefficient models only consider the liquid phase).
In these cases, we need to combine various models together to obtain the 'full' representation.
`CompositeModels` allows users to mix-and-match all of our available models.
In the most general case, five models must be specified:

```julia
struct CompositeModel{𝕍,𝕃,𝕊,𝕃𝕍,𝕃𝕊} <: EoSModel
    components::Vector{String}
    gas::𝕍
    liquid::𝕃
    solid::𝕊
    saturation::𝕃𝕍
    melting::𝕃𝕊
end
```

The simplest case to consider is where we use the ideal gas model to represent the vapour phase, a correlation for the liquid (`RackettLiquid`) and saturation curve (`LeeKeslerSat`):

```julia
julia> model = CompositeModel(["water"])
Composite Model:
 Gas Model: BasicIdeal()
 Liquid Model: RackettLiquid("water")
 Saturation Model: LeeKeslerSat("water")
```

The possibilities with this methodology are truly limitless.
A useful example is in the case of [Solid-liquid equilibrium](./sle_phase_diagrams.md) calculations.

## Group-Contribution Models

Many of the classes of equations of state discussed above also have group-contribution variants.
These methods allow us to assemble species from groups in the cases where pure-component parameters are not available.
Examples include `UNIFAC` and `SAFTgammaMie`.
Let us consider `SAFTgammaMie` first:

```julia
julia> model = SAFTgammaMie([("butane",["CH3"=>2,"CH2"=>2])])
SAFTgammaMie{BasicIdeal, SAFTVRMie{BasicIdeal}} with 1 component:
 "butane": "CH3" => 2, "CH2" => 2
Group Type: SAFTgammaMie
Contains parameters: segment, shapefactor, lambda_a, lambda_r, sigma, epsilon, epsilon_assoc, bondvol

julia> model.params.epsilon
2×2 PairParam{Float64}(["CH3", "CH2"]) with values:
 256.77  350.77
 350.77  473.39

julia> model.groups
GroupParam(:SAFTgammaMie) with 1 component:
 "butane": "CH3" => 2, "CH2" => 2
```

As we can see, we have assembled butane from two methyl and methylene groups.
As such, the parameters within SAFT-$\gamma$ Mie now pertain to the groups, rather than the species.
We also have a new field, `groups`, which provides all the details on the multiplicity of each group.

!!! tip "Available Groups"
    A full list of groups available for each equation of state is available in their respective docs.

### Structured groups

There is an additional class of group contribution model where not only does the group multiplicity have to be specified, but the number of bonds between groups must also be specified.
An example of this is `structSAFTgammaMie`:

```julia
julia> model = sSAFTgammaMie([("butane",["CH3"=>2,"CH2"=>1],[("CH3","CH2")=>2,("CH2","CH2")=>1])])
structSAFTgammaMie{BasicIdeal, SAFTVRMie{BasicIdeal}} with 1 component:
 "butane": "CH3" => 2, "CH2" => 1
Group Type: SAFTgammaMie
Contains parameters: segment, shapefactor, lambda_a, lambda_r, sigma, epsilon, epsilon_assoc, bondvol
```

Note how we must now provide the number of times each group is bonded with one another.
The only other equation of state which requires this is `gcPCSAFT`.

!!! note "External Packages"
    This process can get quite tedious when dealing with larger species.
    As such, we are currently developping a new package to automate this process.
    We hope to register [GCIdentifier.jl](https://github.com/ClapeyronThermo/GCIdentifier.jl) soon.
