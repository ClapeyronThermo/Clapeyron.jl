"""
    salt_stoichiometry(model::ElectrolyteModel,salts)
Obtains the stoichiometry matrix of `salts` made up of ions stored in the `model`. This will also check that the salt is electroneutral and that all ions are involved in the salts.
"""
function salt_stoichiometry(model::ElectrolyteModel,salts)
    iions = model.charge.!=0
    ions = model.components[iions]
    ν = zeros(length(salts),length(ions))
    charges = model.charge[iions]

    for i ∈ 1:length(salts)
        v = salts[i][2]
        for j in 1:length(v)
            name,vj = v[j]
            for k in 1:length(ions)
                if name == ions[k]
                    ν[i,k] = vj
                end
            end
        end
        if dot(@view(ν[i,:]),charges)!==0.
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

function salt_stoichiometry(model::ElectrolyteModel,salts::GroupParam)
    ions = model.components[model.charge.!=0]
    nsalts = length(salts.components)
    ν = zeros(nsalts,length(ions))
    charges = model.charge[model.charge.!=0]
    for i ∈ 1:nsalts
        v = salts.n_groups[i]
        salt_names = salts.groups[i]
        for j in 1:length(v)
            ν[i,salt_names[j].==ions] .= v[j]
        end
        if dot(@view(ν[i,:]),charges)!==0.
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

    nsalts = salts isa GroupParam ? length(salts.components) : length(salts)
    isalts = 1:nsalts
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

function debye_length(V,T,z,ϵ_r,Z)
    s = e_c*e_c/(ϵ_0*ϵ_r*k_B*T)
    I = @sum(z[i]*Z[i]*Z[i])
    κ = Solvers.strong_zero(I) do ii
        sqrt(s*N_A/V)*sqrt(ii)
    end
end

function a_ion(ionmodel, rsp, neutralmodel, V, T, z, neutral_data, ϵ_r)
    return a_ion(ionmodel, V, T, z, ϵ_r)
end

export molality_to_composition
