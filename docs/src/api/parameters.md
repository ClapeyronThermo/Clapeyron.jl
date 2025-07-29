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

## Combining Rules, out-of-place methods

```@docs
Clapeyron.kij_mix
Clapeyron.pair_mix
Clapeyron.mirror_pair
Clapeyron.sigma_LorentzBerthelot
Clapeyron.epsilon_LorentzBerthelot
Clapeyron.epsilon_HudsenMcCoubrey
Clapeyron.epsilon_HudsenMcCoubreysqrt
Clapeyron.lambda_LorentzBerthelot
Clapeyron.lambda_squarewell
```

## Combining Rules, in-place methods

```@docs
Clapeyron.kij_mix!
Clapeyron.pair_mix!
Clapeyron.mirror_pair!
Clapeyron.sigma_LorentzBerthelot!
Clapeyron.epsilon_LorentzBerthelot!
Clapeyron.epsilon_HudsenMcCoubrey!
Clapeyron.epsilon_HudsenMcCoubreysqrt!
Clapeyron.lambda_LorentzBerthelot!
Clapeyron.lambda_squarewell!
```

## Group Combining Rules

```@docs
Clapeyron.group_sum
Clapeyron.group_pairmean
Clapeyron.group_matrix
```

## Model Splitting

```@docs
Clapeyron.split_model
Clapeyron.is_splittable
Clapeyron.index_reduction
```

## Model Exporting

```@docs
Clapeyron.export_model
```

## Model Citing

```@docs
Clapeyron.doi
Clapeyron.cite
Clapeyron.doi2bib
```
