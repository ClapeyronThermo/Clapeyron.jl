# User-Defined Parameters

While Clapeyron has an extensive database for most equations of state, users may need to define their own parameters.
This is possible in Clapeyron whenever constructing your models using the `userlocations` optional argument:

```julia
julia> model = PCSAFT(["your_species"]; userlocations=["your_parameters.csv"])
```

Note that, in cases like cubics where one model can be constructed from multiple other models, you may need to specify the parameters for the submodels separately:

```julia
julia> model = PR(["your_species"]; userlocations=["your_parameters.csv"],
                                    alpha_userlocations=["your_alpha_parameters.csv"])
```

There are two ways which users can specify their own parameters.
The preferred way is to use CSVs.
However, we also allow users to define parameters straight in the REPL.

## Defining parameters using CSVs

By default, parameters in Clapeyron are stored in CSV files.
This was to make them easily readable and modifiable.
You find our [database online](https://github.com/ClapeyronThermo/Clapeyron.jl/tree/master/database) for examples.
Nevertheless, recalling that we allow for three types of parameters (`SingleParam`, `PairParam` and `AssocParam`), we have three types of CSV files:

- Like parameters: These files contain parameters only referring to a single species:

  | Clapeyron Database File      |        |        |
  | ---------------------------- | ------ | ------ |
  | {MODEL NAME} Like Parameters |        |        |
  | species                      | param1 | param2 |
  | water                        | 1234   | 5.678  |

  To ensure that you contain all the appropriate parameters, please check the docs for the respective equation of state.
  The only equation of state that requires additional parameters is SAFT equations where users must also specify the occurrence of each site type (`n_X`) present.

- Unlike parameters: These files contain parameters relating to a pair of species:

  | Clapeyron Database File        |          |       |
  | ------------------------------ | -------- | ----- |
  | {MODEL NAME} Unlike Parameters |          |       |
  | species1                       | species2 | param |
  | water                          | methanol | 0.910 |

  With the exception of activity coefficient models, these parameters are optional.
  Within Clapeyron, if a binary parameter is not specified, a combining rule will be used instead to obtain the unlike parameter.
  Furthermore, if a binary interaction parameter such as $k_{ij}$ (for $\epsilon$ in SAFT or $a$ in cubics) or $l_{ij}$ (for $\sigma$ in SAFT or $b$ in cubics) instead of the raw pair parameters, Clapeyron will use them in the combining rule.
  By default, unlike parameters are treated as symmetric ($ij=ji$).
  In the case of activity coefficient models, unlike parameters must be specified.
  As unlike parameters in activity coefficient models are almost always asymmetric, users must specify the `param_12` and `param_21` parameters in separate rows where the order is determined by which species is `species1` and `species2`.

- Association parameters: These files contain those parameters involved in modelling association interactions:

  | Clapeyron Database File       |       |          |       |       |
  | ----------------------------- | ----- | -------- | ----- | ----- |
  | {MODEL NAME} Assoc Parameters |       |          |       |       |
  | species1                      | site1 | species2 | site2 | param |
  | water                         | H     | water    | e     | 1.234 |
  | water                         | H     | methanol | e     | 5.678 |
  | methanol                      | H     | water    | e     | 9.101 |

  For association parameters, users must specify the species ($ij$) *and* sites ($ab$) involved.
  However, to avoid having to specify a large number of parameters, it is generally assumed that $ij,ab=ji,ba$ which should always be true.
  To allow for a greater degree of flexibility, we do not assume that $ij,ab=ij,ba$ as shown above.
  In Clapeyron, if a set of association parameters isn't specified, it is assumed that those association interactions do not occur.
  The exception to this is if a combining rule is specified within `AssocOptions` (see [relevant docs](./basics_model_construction.md)).

!!! note "Specifying group parameters"
    In group-contribution based approaches, rather than giving the name of the species under the `species` headers, give the group names.
    The rest will work the same way.

If ever you are unsure, please check our database online to see if the structure of your files make sense.

We offer two options to creating CSVs in the REPL:

1. Using `ParamTable`: you can specify `:single`, `:pair` and `:assoc` tables which just need to match the column layout in the original tables:

```julia
data = (species = ["water"], Mw = [18.0])
file = ParamTable(:single,data,name="water_mw")
model = PR(["water"],user_locations = [file])
```

1. Using raw CSVs: you can literally write the CSV out as a string within the REPL and just substitute this in:

```julia
csv_data = """Clapeyron Database File,
       my water like parameters
       species,Mw
       water,18.0
       """
model = PR(["water"],user_locations = [csv_data])
```

We prefer that users use the CSV approach as they are much easier to visualise and modify.
However, it is a more tedious approach.
As such, we provide a more direct way described below.

## Defining parameters within the REPL

Still using the `userlocations` argument, we can specify the parameters used directly within the REPL.
For example, in the case of a cubic:

```julia
julia> model = PR(["water","methanol"]; userlocations = (;
                    Tc = [647.13,512.64],
                    Pc = [2.19e7,8.14e6],
                    Mw = [18.0,32.0],
                    k  = [ 0.00 -0.01;
                          -0.01  0.00]),
                                  alpha_userlocations = (;
        acentricfactor = [0.343,0.566]))
PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule} with 2 components:
 "water"
 "methanol"
Contains parameters: a, b, Tc, Pc, Mw
```

Single-component parametes are specified as vectors and pair parameters are specified as matrices.
Note that you do still need to specify the submodel parameters separately (in the case of the alpha function above).

As an additional example, let us now consider a similar case in a SAFT-type equation of state:

```julia
model = SAFTVRMie(["water","methanol"]; userlocations=(;
       Mw            = [18.01,32.04],
       segment       = [1.0,1.67034],
       epsilon       = [266.68 278.45;
                        278.45 307.69],
       sigma         = [3.0063,3.2462],
       lambda_a      = [6.,6.],
       lambda_r      = [17.02,7.6134],
       n_H           = [2,1],
       n_e           = [2,2],
       epsilon_assoc = Dict((("water","e"),("water","H")) => 1985.4,
                            (("methanol","e"),("methanol","H")) => 2062.1,
                            (("methanol","e"),("water","H")) => 1993.5,
                            (("water","e"),("methanol","H")) => 1993.5),
       bondvol       = Dict((("water","e"),("water","H")) => 1.0169e-28,
                            (("methanol","e"),("methanol","H")) => 1.0657e-28,
                            (("methanol","e"),("water","H")) => 1.0411e-28,
                            (("water","e"),("methanol","H")) => 1.0411e-28)))
```

As we can see, with certain equations of state, this method of specifying parameters can become unwieldy.
This is why we recommend using CSVs instead.
Also note that, if you only specify pure component parameters for a parameter that should include pair parameters (such as `sigma` above), combining rules will be used to obtain the remaining parameters.
If you specify all of the pair parameters, they will be used instead (in the case of `epsilon`).
For the cubic equation of state earlier, we specified `k`, which will be used to obtain the unlike `a` parameters.
Further, if we had not specified the cross-association parameters between methanol and water, unless we specified otherwise in `AssocOptions`, these interactions would not be included.
