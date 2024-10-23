```@meta
CurrentModule = Clapeyron
```

## Index

```@index
Pages = ["dev.md"]
```

# Developer guide

## Deactivating precompilation

From Julia 1.9 onwards, the precompilation possibilities of Julia are vastly improved, so we opt to precompile some frequently used methods. this can be disabled by calling `precompile_clapeyron!(false)`

```@docs
Clapeyron.precompile_clapeyron!
```
