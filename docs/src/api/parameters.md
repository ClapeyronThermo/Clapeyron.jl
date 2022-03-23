```@meta
CurrentModule = Clapeyron
```

## Contents

```@contents
Pages = ["parameters.md"]
```

## Index

```@index
Pages = ["parameters.md"]
```

## Parsing Parameters from Files
```@docs
Clapeyron.ParamOptions
Clapeyron.getparams
```

## Creating Files from Parameters
```@docs
Clapeyron.ParamTable
Clapeyron.cleartemp!
```

## Parameter types
```@docs
Clapeyron.SingleParam
Clapeyron.PairParam
Clapeyron.AssocParam
Clapeyron.GroupParam
Clapeyron.SiteParam
Clapeyron.AssocOptions
```

## Combining Rules
```
Clapeyron.kij_mix
Clapeyron.pair_mix
Clapeyron.sigma_LorentzBerthelot
Clapeyron.epsilon_LorentzBerthelot
Clapeyron.epsilon_HudsenMcCoubrey
Clapeyron.lambda_LorentzBerthelot
Clapeyron.lambda_squarewell
Clapeyron.group_sum
```

## Model Splitting
```@docs
Clapeyron.split_model
Clapeyron.is_splittable
Clapeyron.index_reduction
```
