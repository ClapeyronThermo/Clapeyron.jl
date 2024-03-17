struct EoSVirial2{I,T} <: SecondVirialModel
    idealmodel::I
    model::T
end

"""
    EoSVirial2 <: SecondVirialModel
    EoSVirial2(model;idealmodel = idealmodel(model))

## Input models
- `model`: Model providing second virial coefficient is obtained
- `idealmodel`: Ideal Model

## Description

Virial model, that just calls the second virial coefficient of the underlying model.
```
B(T,z) = B(model,T,z)
```
"""
EoSVirial2

function EoSVirial2(model;idealmodel = idealmodel(model))
    return EoSVirial2{typeof(idealmodel),typeof(model)}(idealmodel,model)
end

second_virial_coefficient_impl(model::EoSVirial2,T,z=SA[1.0]) = second_virial_coefficient_impl(model.model,T,z)

export EoSVirial2