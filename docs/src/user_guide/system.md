# Developing OpenSAFT
There are three components to the implementation of an equation of state (EoS):

1. importing parameters, and getting them to a form where they can easily be accessed
2. writing the actual EoS, and
3. generating physical properties out of the EoS.

## 1. Importing parameters
OpenSAFT is designed such that raw data files can systematically extracted so that the developer will be able to focus on writing the EoS. The database extraction is handled by the functions in ```utils/extractdatabase.jl```.  A full documentation of how this is set up will be written in the near future, but the details of the implementation are not necessary for the EoS developer.

Currently, the raw data files are formatted as CSVs with a semicolon delimiter (techinically DSVs), with the headers on the third row. We currently support three parameter data files for each EoS:

1. like parameters, for parameters that involve a single component
2. unlike parameters, for parameters that involve two components, and
3. assoc parameters, for parameters involved in the association term.

We may support more data files in the future, such as pure heat capacity polynomial coefficients for a more sophisticated ideal term.

Column data in the relevant rows are extracted using the function

```julia
retrieve_params(components, selected method;
    customdatabase_like = "None", variant_like = "None", redirect_like = "None",
    customdatabase_unlike = "None", variant_unlike = "None", redirect_unlike = "None",
    customdatabase_assoc = "None", variant_assoc = "None", redirect_assoc = "None")
```

By default, it looks for database files for the relevant method in the ```database``` directory of the project root directory, but you may choose to redirect the search to a different directory using optional parameters. This function returns an Array of dictionaries, one for each database file.

Once the column data are extracted, parameter data can be filtered using

```julia
function filterparams(raw_params, like_params;
    unlike_params = [], assoc_params = [])
```

The ```like_params```, ```unlike_params```, and ```assoc_params``` is a list of column headers to extract to their individidual dictionaries.

The keys for the main species are Sets, because they have the property where ```Set([i,j]) == Set([j,i])```, and ```Set([i,i]) == Set([i])```, which is perfect for the representation of symmetric interactions that we observe in thermodynamics, and it simplifies loops in the equations.

For interaction parameters between species ```i``` and ```j```, they are accessed by the key ```Set([i,j])```, while interaction parameters between the same species ```i``` are accessed by the key ```Set([i])```.

For the association parameters, the keys are ```Set([(i,a),(j,b)])```, where ```i``` and ```j``` are the main species, and ```a``` and ```b``` are the asociation species on species ```i``` and ```j``` respectively.

## 2. Writing the EoS
OpenSAFT makes use of Julia's multiple dispatch to reduce function duplication for variant SAFT models of the same family. The Structs are defined in ```models/model_structs.jl```, where the hierachy is as follows:

```julia
{SAFT_model} <: {SAFT_family} <: SAFT <: EoS
```

In the future, if we ever wish to provide support for non-SAFT EoS (cubics, activity coefficient based, etc), we may add a new abstract type under EoS.

Parameters that are used by the model are defined in ```models/param_structs.jl```.

The user will call the ```system(components, method; kwargs)``` function in ```models/system.jl``` to define their system. This function extracts the raw data using function ```retrieveparams```, and these extracted parameters are interpreted in the function ```create_{model}Params``` in the file ```models/import_params.jl```.

The EoS for SAFT is expressed in terms of the Helmholtz free energy, and this is defined in ```models/eos/eos.jl```. This calls the relevant functions to generate, for example ```a_ideal``` and ```a_res```. For SAFT, these functions will reside in the ```models/eos/SAFT``` directory.

## 3. Generating physical properties out of the EoS
All methods are stored in the ```methods``` directory. Currently, we have support for getting pure properties for SAFT. Additional methods to support mixtures, etc, will be written here. We may make use of Julia's multiple dispatch to create model-specific methods if that is necessary.

# Solvers
There may be a few solvers that are typically used in thermodynamics but are not currently readily available in the Julia community, so they will be set up ande put together here in the ```solvers``` directory in a module called Solvers.
