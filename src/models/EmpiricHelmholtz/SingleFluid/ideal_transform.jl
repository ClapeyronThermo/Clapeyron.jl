"""
    idealmodel_to_json_data(model::EoSModel;Tr = 1.0,T0 = 298.15, Vr = 1.0)

Transforms an `model::IdealModel` into a vector of dictionaries containing valid ideal multiparameter helmholtz terms.
`Tr` is the reducing temperature, `T0` is the reference temperature, `Vr` is the reducing volume.
## Example 
```
julia> id = BasicIdeal(["water"])
BasicIdeal(Clapeyron.BasicIdealParam)

julia> Clapeyron.idealmodel_to_json_data(id)
1-element Vector{Dict{Symbol, Any}}:
 Dict(:T0 => 298.15, :type => "IdealGasHelmholtzCP0Constant", :cp_over_R => 2.5, :Tc => 1.0)

```
"""
function idealmodel_to_json_data(model;Tr = 1.0,T0 = 298.15,Vr = 1.0)
    return idealmodel_to_json_data(model,Tr,T0,Vr)
end

function idealmodel_to_json_data(model::BasicIdealModel,Tr,T0,Vr)
    [
        Dict(:type => "IdealGasHelmholtzCP0Constant",
        :cp_over_R => 1.5,
        :T0 => T0,
        :Tc => Tr
        ),
    ]
end

function idealmodel_to_json_data(model::ReidIdealModel,Tr,T0,Vr)
    single_component_check(idealmodel_to_json_data,model)
    coeffs = model.params.coeffs[1]
    [
        Dict(:type => "IdealGasHelmholtzCP0PolyT",
        :c => [coeffs...],
        :t => [0,1,2,3],
        :Tc => Tr,
        :T0 => T0,
        )
    ]
end

function idealmodel_to_json_data(model::JobackIdealModel,Tr,T0,Vr)
    single_component_check(idealmodel_to_json_data,model)
    return idealmodel_to_json_data(ReidIdeal(model),Tr,T0,Vr)
end

function idealmodel_to_json_data(model::MonomerIdealModel,Tr,T0,Vr)
    single_component_check(idealmodel_to_json_data,model)
    Mwᵢ = model.params.Mw[1]*0.001
    Λᵢ = h/√(k_B*Mwᵢ/N_A) # * T^(-1/2)
    kᵢ = N_A*Λᵢ^3 #T^(-3/2)
    # monomer: a = ∑ xi * [log(xi*ki*T^-1.5/v)]
    # ∑ xi * [log(xi) +  1.5*log(ki*T/v)]
    # ∑ xi * [log(xi) +  a0i(v,T)]
    #a0i(v,T) = log(ki) - log(v) + 1.5*log(Tinv)
    #a0i(v,T) = log(ki) - log(v) + log(vr) - log(vr) + 1.5*log(Tinv) + 1.5*log(Tr) - 1.5*log(Tr)
    #a0i(v,T) = log(ki) + log(vr/v) - log(vr)  - 1.5*log(Tr) + 1.5*log(Tr/Tinv)
    #a0i(v,T) = log(vr/v)  + log(ki)- log(vr) - 1.5*log(Tr) + 1.5*log(Tr/Tinv)
    #a1 = log(ki) - log(vr) - 1.5*log(Tr)
    #a2 = 1.5
    a1 = log(kᵢ) - log(Vr) - 1.5*log(Tr)
    [
        Dict(:type => "IdealGasHelmholtzLead",
            :a1 => a1,
            :a2 => 0.0,
        )
        Dict(:type => "IdealGasHelmholtzLogTau",
            :a => 1.5,
        )
    ]
end






