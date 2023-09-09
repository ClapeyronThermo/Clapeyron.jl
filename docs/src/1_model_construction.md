# 1. Model Construction
All functions in Clapeyron revolve around an object we refer to as the `model`. These models are intended to hold all of the information required to model a given system within a specific equation of state. The absolute simplest of these is the `BasicIdeal` model, which models species as an ideal gas with only translational modes of motion. It can be constructed simply using:
```julia
julia> model = BasicIdeal(["water"])
BasicIdeal()
```
This `model` is unique as it is the only model object that does not hold any information as the ideal gas is a universal model that does not depend on the chemical identity of the species. Nevertheless, one can use this model to obtain properties such as the volume and isobaric heat capacity at a given temperature and pressure:
```julia
julia> volume(model,1e5,298.15)
0.02478957029602388

julia> isobaric_heat_capacity(model,1e5,298.15)
20.786156545383097
```
At this stage, we point out that all values within Clapeyron are in SI units, unless otherwise specified.

From here, we will now consider how `model`s are constructed in different equations of state and how users can implement their own parameters.

## 1.1 Ideal gas models

## 1.2 Cubic equations of state

## 1.3 SAFT equations of state

## 1.4 Empirical equations of state

## 1.5 Activity coefficient models

## 1.6 Composite models

## 1.A List of available groups
### 1.A.1 Joback Method
### 1.A.2 Walker Ideal Gas Method
### 1.A.3 SAFT-$\gamma$ Mie
### 1.A.4 gc-PC-SAFT
### 1.A.5 UNIFAC
#### 1.A.5.1 (mod) UNIFAC
#### 1.A.5.2 og-UNIFAC
