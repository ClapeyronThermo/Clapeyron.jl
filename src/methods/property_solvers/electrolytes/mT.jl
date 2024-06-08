"""
    mean_ionic_activity_coefficient_sat(model::ESElectrolyteModel,salts,T,m,zsolvent=[1.])
Calculate the mean ionic activity coefficient of selection of salts at the saturation point at a certain temperature and molality. These are defined as:
```
γ± = φ±/φ±₀ * ∑zsolv/∑z
```
Example:
```julia    
model = ePCSAFT(["water"],["sodium","chloride"])

salts = [("sodium chloride",("sodium"=>1,"chloride"=>1))]

T = 298.15
m = [1.0]

γ± = mean_ionic_activity_coefficient_sat(model,salts,T,m)
```
If multiple solvents are present, the composition of the solvent can be specified with the `zsolvent` keyword argument.
"""
function mean_ionic_activity_coefficient_sat(model::ESElectrolyteModel,salts,T,m,zsolvent=[1.])
    isolvent = model.icomponents[model.charge.==0]
    iions = model.icomponents[model.charge.!=0]
    ions = model.components[model.charge.!=0]

    method = FugBubblePressure(nonvolatiles=ions)

    ν = salt_stoichiometry(model,salts)

    z0 = molality_to_composition(model,salts,ones(length(m)).*1e-20,zsolvent)
    z = molality_to_composition(model,salts,m,zsolvent)

    (p,vl,vv,y) = bubble_pressure(model,T,z,method)
    φ = fugacity_coefficient(model,p,T,z;phase=:l)[iions]
    φ0 = fugacity_coefficient(model,p,T,z0;phase=:l)[iions]

    γim = φ./φ0.*sum(z[isolvent])/sum(z)
    γsm = (prod(γim'.^ν,dims=2)).^(1 ./sum(ν,dims=2))
    return γsm
end

"""
    mean_ionic_activity_coefficient(model::ESElectrolyteModel,salts,p,T,m,zsolvent=[1.])
Calculate the mean ionic activity coefficient of selection of salts at a given pressure, temperature and molality. These are defined as:
```
γ± = φ±/φ±₀ * ∑zsolv/∑z
```
Example:
```julia    
model = ePCSAFT(["water"],["sodium","chloride"])

salts = [("sodium chloride",("sodium"=>1,"chloride"=>1))]

p = 1e5
T = 298.15
m = [1.0]

γ± = mean_ionic_activity_coefficient(model,salts,p,T,m)
```
If multiple solvents are present, the composition of the solvent can be specified with the `zsolvent` keyword argument.
"""
function mean_ionic_activity_coefficient(model::ESElectrolyteModel,salts,p,T,m,zsolvent=[1.])
    isolvent = model.icomponents[model.charge.==0]
    iions = model.icomponents[model.charge.!=0]
    ions = model.components[model.charge.!=0]

    ν = salt_stoichiometry(model,salts)

    z0 = molality_to_composition(model,salts,ones(length(m)).*1e-20,zsolvent)
    z = molality_to_composition(model,salts,m,zsolvent)

    φ0 = fugacity_coefficient(model,p,T,z0)[iions]

    φ = fugacity_coefficient(model,p,T,z)[iions]

    γim = φ./φ0.*sum(z[isolvent])/sum(z)
    γsm = (prod(γim'.^ν,dims=2)).^(1 ./sum(ν,dims=2))
    return γsm
end


"""
    osmotic_coefficient_sat(model::ESElectrolyteModel,salts,T,m,zsolvent=[1.])
Calculate the osmotic coefficient of selection of solvents at the saturation point at a certain temperature and molality. These are defined as:
```
ϕ = -1/(∑νi*mi*Mw)*log(asolv)
```
Example:
```julia    
model = ePCSAFT(["water"],["sodium","chloride"])

salts = [("sodium chloride",("sodium"=>1,"chloride"=>1))]

T = 298.15
m = [1.0]

ϕ = osmotic_coefficient(model,salts,T,m)
```
If multiple solvents are present, the composition of the solvent can be specified with the `zsolvent` keyword argument.
"""
function osmotic_coefficient_sat(model::ESElectrolyteModel,salts,T,m,zsolvent=[1.])
    isolvent = model.icomponents[model.charge.==0]
    iions = model.icomponents[model.charge.!=0]
    ions = model.components[model.charge.!=0]

    method = FugBubblePressure(nonvolatiles=ions)

    ν = salt_stoichiometry(model,salts)

    z0 = molality_to_composition(model,salts,ones(length(m)).*1e-20,zsolvent)
    z = molality_to_composition(model,salts,m,zsolvent)

    (p,vl,vv,y) = bubble_pressure(model,T,z,method)
    φ = fugacity_coefficient(model,p,T,z;phase=:l)[isolvent]
    φ0 = fugacity_coefficient(model,p,T,z0;phase=:l)[isolvent]

    asolv = φ./φ0.*z[isolvent]/sum(z)
    Mw = mw(model.neutralmodel)[isolvent].*1e-3
    # println()
    return -1 ./(sum(ν.*m).*Mw).*log.(asolv)
end

"""
    osmotic_coefficient(model::ESElectrolyteModel,salts,p,T,m,zsolvent=[1.])
Calculate the osmotic coefficient of selection of solvents at a given pressure, temperature and molality. These are defined as:
```
ϕ = -1/(∑νi*mi*Mw)*log(asolv)
```
Example:
```julia    
model = ePCSAFT(["water"],["sodium","chloride"])

salts = [("sodium chloride",("sodium"=>1,"chloride"=>1))]

p = 1e5
T = 298.15
m = [1.0]

ϕ = osmotic_coefficient(model,salts,p,T,m)
```
If multiple solvents are present, the composition of the solvent can be specified with the `zsolvent` keyword argument.
"""
function osmotic_coefficient(model::ESElectrolyteModel,salts,p,T,m,zsolvent=[1.])
    isolvent = model.icomponents[model.charge.==0]

    ν = salt_stoichiometry(model,salts)

    z0 = molality_to_composition(model,salts,ones(length(m)).*1e-20,zsolvent)
    z = molality_to_composition(model,salts,m,zsolvent)

    φ0 = fugacity_coefficient(model,p,T,z0;phase=:l)[isolvent]

    φ = fugacity_coefficient(model,p,T,z;phase=:l)[isolvent]

    asolv = φ./φ0.*z[isolvent]/sum(z)
    Mw = mw(model.neutralmodel)[isolvent].*1e-3
    # println()
    return -1 ./(sum(ν.*m).*Mw).*log.(asolv)
end

export mean_ionic_activity_coefficient, osmotic_coefficient, mean_ionic_activity_coefficient_sat, osmotic_coefficient_sat