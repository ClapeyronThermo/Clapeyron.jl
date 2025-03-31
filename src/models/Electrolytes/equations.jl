"""
    salt_stoichiometry(model::ElectrolyteModel,salts)
Obtain the stoichiometry matrix of `salts` made up of ions stored in the `model`. This will also check that the salt is electroneutral and that all ions are involved in the salts.
"""
function salt_stoichiometry(model::ElectrolyteModel,salts)
    ions = model.components[model.charge.!=0]
    ν = zeros(length(salts),length(ions))
    for i ∈ 1:length(salts)
        v = salts[i][2]
        for j in 1:length(v)
            ν[i,v[j][1].==ions] .= v[j][2]
        end
        if sum(ν[i,:].*model.charge[model.charge.!=0])!==0.
            throw(ArgumentError("The salt $i is not electroneutral"))
        end
    end
    for νi in eachcol(ν)
        if iszero(sum(νi))
            throw(ArgumentError("Not all ions are involved in the salts"))
        end
    end
    return ν
end

"""
    molality_to_composition(model::ElectrolyteModel,salts,m,zsolv=[1.])

Convert molality (mol/kg) to composition for a given model, salts, molality, and solvent composition.
"""
function molality_to_composition(model::ElectrolyteModel,salts,m,zsolv=SA[1.],ν = salt_stoichiometry(model,salts))
    nc = length(model)
    Z = model.charge
    nions = count(!iszero,Z)
    nneutral = nc - nions
    Mw = mw(model.neutralmodel).*1e-3
    isalts = 1:length(salts)
    iions = 1:nions
    ineutral = 1:nneutral
    ∑mν = sum(m[k]*sum(ν[k,i] for i ∈ iions) for k ∈ isalts)
    x_solv = zsolv ./ (1+sum(zsolv[j]*Mw[j] for j in ineutral)*∑mν)
    x_ions = [sum(m[k]*ν[k,l] for k ∈ isalts) / (1/sum(zsolv[j]*Mw[j] for j in ineutral)+∑mν) for l ∈ iions]
    return vcat(x_solv,x_ions)
end

"""
    @iions()

This macro is a non-allocating equivalent to the following code:

```julia
(1:length(model))[model.params.charges.values .!= 0]
```

`@iions` is an iterator that goes through all charged components in an electrolyte model.
"""
macro iions()
    quote
        Iterators.filter(!iszero ∘ Base.Fix1(Base.getindex,Z), 1:length(Z))
    end |> esc
end

"""
    @ineutral()

This macro is a non-allocating equivalent to the following code:

```julia
(1:length(model))[model.params.charges.values .== 0]
```

`@iions` is an iterator that goes through all non charged components in an electrolyte model.
"""
macro ineutral()
    quote
        Iterators.filter(iszero ∘ Base.Fix1(Base.getindex,Z), 1:length(Z))
    end |> esc
end

function debye_length(V,T,z,ϵ_r,Z,∑z = sum(z))
    ρ = N_A*∑z/V
    s = e_c*e_c/(4π*ϵ_0*ϵ_r*k_B*T)
    κ = sqrt(4π*s*ρ*@sum(z[i]*Z[i]*Z[i])/∑z)
end

function a_ion(ionmodel, rsp, neutralmodel, V, T, z, neutral_data, ϵ_r)
    return a_ion(ionmodel, V, T, z, ϵ_r)
end
export molality_to_composition