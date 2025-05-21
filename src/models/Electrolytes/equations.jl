"""
    salt_stoichiometry(model::ElectrolyteModel,salts)
Obtain the stoichiometry matrix of `salts` made up of ions stored in the `model`. This will also check that the salt is electroneutral and that all ions are involved in the salts.
"""
function salt_stoichiometry(model::ElectrolyteModel,salts)
    ions = model.components[model.charge.!=0]
    ν = zeros(length(salts),length(ions))
    for i in 1:length(salts)
        v = salts[i][2]
        for j in 1:length(v)
            ν[i,v[j][1].==ions] .= v[j][2]
        end
        if sum(ν[i,:].*model.charge[model.charge.!=0])!==0.
            throw(ArgumentError("The salt $i is not electroneutral"))
        end
    end
    if any(sum(ν;dims=1).==0)
        throw(ArgumentError("Not all ions are involved in the salts"))
    end

    return ν
end

"""
    molality_to_composition(model::ElectrolyteModel,salts,m,zsolv=[1.])

Convert molality (mol/kg) to composition for a given model, salts, molality, and solvent composition.
"""
function molality_to_composition(model::ElectrolyteModel,salts,m,zsolv=SA[1.],ν = salt_stoichiometry(model,salts))
    ions = model.components[model.charge.!=0]
    neutral = model.components[model.charge.==0]
    Mw = mw(model.neutralmodel).*1e-3
    isalts = 1:length(salts)
    iions = 1:length(ions)
    ineutral = 1:length(neutral)
    ∑mν = sum(m[k]*sum(ν[k,i] for i ∈ iions) for k ∈ isalts)
    x_solv = zsolv ./ (1+sum(zsolv[j]*Mw[j] for j in ineutral)*∑mν)
    x_ions = [sum(m[k]*ν[k,l] for k ∈ isalts) / (1/sum(zsolv[j]*Mw[j] for j in ineutral)+∑mν) for l ∈ iions]
    return vcat(x_solv,x_ions)
end

export molality_to_composition