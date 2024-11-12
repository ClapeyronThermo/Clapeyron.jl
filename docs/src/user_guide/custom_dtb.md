## Definitions

Within Clapeyron, we use CSV files to store our parameters.
There are four types of database files for different types of parameters:

- Like parameters: These CSV files have the structure:

  | Clapeyron Database File      |        |        |
  | ---------------------------- | ------ | ------ |
  | {MODEL NAME} Like Parameters |        |        |
  | species                      | param1 | param2 |
  | water                        | 1234   | 5.678  |

  These are used for parameters which only refer to a single species (such as the critical temperature, molar mass, like segment size, number of segments).
  These types of files are also used for the ideal model databases (since all of these are just for like species).

- Unlike parameters: These CSV files have the structure:

  | Clapeyron Database File        |          |       |
  | ------------------------------ | -------- | ----- |
  | {MODEL NAME} Unlike Parameters |          |       |
  | species1                       | species2 | param |
  | water                          | methanol | 0.910 |

  These are used for parameters which refer to a pair of different species (such as the unlike interaction parameter $k_{ij}$).

- Association parameters: These CSV files have the structure:

  | Clapeyron Database File       |       |          |       |       |
  | ----------------------------- | ----- | -------- | ----- | ----- |
  | {MODEL NAME} Assoc Parameters |       |          |       |       |
  | species1                      | site1 | species2 | site2 | param |
  | water                         | H     | water    | e     | 1.234 |
  | water                         | H     | methanol | e     | 5.678 |
  | methanol                      | H     | water    | e     | 5.678 |

  These are used for parameters which refer to a pair of species and sites (such as the association potential depth, `epsilon_assoc`, and bonding volume, `bondvol`).
  Note that this can be for associations between the same species and different sites, or different species and different sites (as shown above).

Note that it is extremely important that the cell A2 has the word 'Like', 'Unlike' or 'Assoc' in it so that Clapeyron can identify the type of parameters in it.

Feel free to check these out in the package to see some better examples!

!!! note "Association Values are asymmetric!"
    There are cases where the association values are asymmetric (for example, water and acetonitrile), so if you are adding cross-association values for a component-site pair, remember to add the corresponding swapped value, like in the water–methanol example above.

## Using your own parameters

If you have CSV files formatted as above with your own parameters, and you want to implement these into one of the existing equations of state in Clapeyron, all that is needed is to provide the path to those files in the definition of your model (note that ideal term related parameters are specified separately):

```julia
model1 = PR(["your_species_1","your_species_2"];userlocations=["path/to/your/database/"], ideal_userlocations=["path/to/your/ideal_database"])
model2 = PCSAFT(["your_species_1","your_species_2"];userlocations=["dtb_like","dtb_unlike","dtb_assoc"],ideal_userlocations=["dtb_ideal"])
```

The rest works exactly as it normally would! We recommend reading the background documentation for the various models, as well as the [`getparams`](@ref) docs, to ensure the units of the parameters you provide are correct and how those parameters are parsed into each model.

You can create those parameters without leaving the Julia REPL, by using [`Clapeyron.ParamTable`](@ref).
This function will create a temporary location on where a CSV containing the table is created:

```julia
data = (species = ["water"],Mw = [18.0])
file = ParamTable(:single,data,name="water_mw")
model = PR(["water"],user_locations = [file])
model.params.Mw[1]  # 18.0
```

You can also write a CSV as a string an pass that directly:

```julia
# Any string that starts with `Clapeyron Database File` will be parsed as a CSV file directly.
csv_data = """Clapeyron Database File,
       my water like parameters
       species,Mw
       water,19.0
       """

model = PR(["water"],user_locations = [csv_data])
model.params.Mw[1]  # 19.0
```
