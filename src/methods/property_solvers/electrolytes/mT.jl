function __electrolyte_fugacities(model,salts,p,T,m,zsolv = SA[1.0];sat = false)
    icomponents = 1:length(model)
    isolvent = icomponents[model.charge.==0]
    iions = icomponents[model.charge.!=0]
    
    ν = salt_stoichiometry(model,salts)
    z0 = molality_to_composition(model,salts,ones(length(m)).*1e-20,zsolv,ν)
    z = molality_to_composition(model,salts,m,zsolv,ν)
    if sat
        (px,vl,vv,y) = bubble_pressure(model,T,z)
        v = vl
    else
        px = p
        v = volume(model,px,T,z;phase=:l)
    end

    RT = Rgas(model)*T
    μ = VT_chemical_potential_res(model,v,T,z)
    Z = px*v/RT/sum(z)
    v0 = volume(model,px,T,z0; phase=:l)
    μ0 = VT_chemical_potential_res(model,v0,T,z0)
    Z0 = px*v0/RT/sum(z0)
    γ = @. exp((μ-μ0)/RT)*Z0/Z
    return (z,z0),γ,(isolvent,iions),ν
end


"""
    mean_ionic_activity_coefficient_sat(model::ESElectrolyteModel,T,m,zsolv = SA[1.0])
    mean_ionic_activity_coefficient_sat(model::ESElectrolyteModel,salts,T,m,zsolv = SA[1.0])

Calculates the mean ionic activity coefficient of selection of salts at the saturation point at a certain temperature `T` and molality `m`. 
These are defined as:

```
γ± = φ±/φ±₀ * ∑zsolv/∑z
```

## Example:

```julia    
model = ePCSAFT(["water"],["sodium","chloride"])

T = 298.15
m = [1.0]

#supposing all combinations of binary salts
γ± = mean_ionic_activity_coefficient_sat(model,T,m)

#specify your own salt composition
salts = [("sodium chloride",("sodium"=>1,"chloride"=>1))]
γ± = mean_ionic_activity_coefficient_sat(model,salts,T,m)
```
If multiple solvents are present, the composition of the solvent can be specified with the `zsolv` keyword argument.
"""
function mean_ionic_activity_coefficient_sat(model::ESElectrolyteModel,salts::Union{AbstractVector,GroupParam},T::Number,m::Number,zsolv=SA[1.])
    return mean_ionic_activity_coefficient(model,salts,zero(T),T,m,zsolv,true)
end

function mean_ionic_activity_coefficient_sat(model::ESElectrolyteModel,T::Number,m::Number,zsolv = SA[1.0])
    return mean_ionic_activity_coefficient_sat(model,auto_binary_salts(model),T,m,zsolv)
end

"""
    mean_ionic_activity_coefficient(model::ESElectrolyteModel,p,T,m,zsolv = SA[1.0])
    mean_ionic_activity_coefficient(model::ESElectrolyteModel,salts,p,T,m,zsolv = SA[1.0])


Calculates the mean ionic activity coefficient of selection of salts at a given pressure `p`, temperature `T` and molality `m`. These are defined as:
```
γ± = φ±/φ±₀ * ∑zsolv/∑z
```
Example:
```julia    
model = ePCSAFT(["water"],["sodium","chloride"])

p = 1e5
T = 298.15
m = [1.0]

#supposing all combinations of binary salts
γ± = mean_ionic_activity_coefficient(model,p,T,m)

#specify your own salt composition
salts = [("sodium chloride",("sodium"=>1,"chloride"=>1))]
γ± = mean_ionic_activity_coefficient(model,salts,p,T,m)
```
If multiple solvents are present, the composition of the solvent can be specified with the `zsolv` keyword argument.
"""
function mean_ionic_activity_coefficient(model::ESElectrolyteModel,salts::Union{AbstractVector,GroupParam},p,T,m,zsolv=SA[1.],sat = false)
    (z,z0),(_γ),(isolvent,iions),ν = __electrolyte_fugacities(model,salts,p,T,m,zsolv,sat = sat)
    γ = _γ[iions]
    γim = γ .* sum(z[isolvent])/sum(z)
    γsm = (prod(γim'.^ν,dims=2)).^(1 ./sum(ν,dims=2))
    return γsm
end

function mean_ionic_activity_coefficient(model::ESElectrolyteModel,p::Number,T,m,zsolv=SA[1.],sat = false)
    return mean_ionic_activity_coefficient(model,auto_binary_salts(model),p,T,m,zsolv,sat)
end

"""
    osmotic_coefficient_sat(model::ESElectrolyteModel,T,m,zsolv = SA[1.0])
    osmotic_coefficient_sat(model::ESElectrolyteModel,salts,T,m,zsolv = SA[1.0])


Calculates the osmotic coefficient of selection of solvents at the saturation point at a certain temperature `T` and molality `m`. These are defined as:
```
ϕ = -1/(∑νi*mi*Mw)*log(asolv)
```
Example:
```julia    
model = ePCSAFT(["water"],["sodium","chloride"])

T = 298.15
m = [1.0]

#supposing all combinations of binary salts
ϕ = osmotic_coefficient_sat(model,T,m)

#specify your own salt composition
salts = [("sodium chloride",("sodium"=>1,"chloride"=>1))]
ϕ = osmotic_coefficient_sat(model,salts,T,m)
```
If multiple solvents are present, the composition of the solvent can be specified with the `zsolv` keyword argument.
"""
function osmotic_coefficient_sat(model::ESElectrolyteModel,salts::Union{AbstractVector,GroupParam},T,m,zsolv = SA[1.0])
   return osmotic_coefficient(model,salts,zero(T),T,m,zsolv,true)
end

osmotic_coefficient_sat(model,T,m,zsolv = SA[1.0]) = osmotic_coefficient_sat(model,auto_binary_salts(model),T,m,zsolv)

"""
    osmotic_coefficient(model::ESElectrolyteModel,salts,p,T,m,zsolv = SA[1.0])
Calculates the osmotic coefficient of selection of solvents at a given pressure `p`, temperature `T` and molality `m`. These are defined as:
```
ϕ = -1/(∑νi*mi*Mw)*log(asolv)
```
Example:
```julia    
model = ePCSAFT(["water"],["sodium","chloride"])

p = 1e5
T = 298.15
m = [1.0]

#supposing all combinations of binary salts
ϕ = osmotic_coefficient(model,p,T,m)

#specify your own salt composition
salts = [("sodium chloride",("sodium"=>1,"chloride"=>1))]
ϕ = osmotic_coefficient(model,salts,p,T,m)
```
If multiple solvents are present, the composition of the solvent can be specified with the `zsolv` keyword argument.
"""
function osmotic_coefficient(model::ESElectrolyteModel,salts::Union{AbstractVector,GroupParam},p,T,m,zsolv=SA[1.0],sat = false)
    (z,z0),(_γ),(isolvent,iions),ν = __electrolyte_fugacities(model,salts,p,T,m,zsolv,sat = sat)
    γ = _γ[isolvent]
    asolv = γ.*z[isolvent]/sum(z)
    Mw = mw(model.neutralmodel)[isolvent].*1e-3
    return -1 ./(sum(ν.*m).*Mw).*log.(asolv)
end

osmotic_coefficient(model::ESElectrolyteModel,p,T,m,zsolv=SA[1.0],sat = false) = osmotic_coefficient(model,auto_binary_salts(model),p,T,m,zsolv,sat)

export mean_ionic_activity_coefficient, osmotic_coefficient, mean_ionic_activity_coefficient_sat, osmotic_coefficient_sat