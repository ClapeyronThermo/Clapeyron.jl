"""
    salt_stoichiometry(model::ElectrolyteModel)
    salt_stoichiometry(model::ElectrolyteModel,salts)

Obtains the stoichiometry matrix of `salts` made up of ions stored in the `model`. 
This will also check that the salt is electroneutral and that all ions are involved in the salts.
If no `salts` argument is specified, it will be created via `Clapeyron.auto_binary_salts(model)`

"""
function salt_stoichiometry(model::ElectrolyteModel,salts = auto_binary_salts(model))
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
    molality_to_composition(model::ElectrolyteModel,salts,m,zsolv = SA[1.0])
    molality_to_composition(model::ElectrolyteModel,m,zsolv = SA[1.0])

Convert molality (mol/kg) to composition for a given model, salts, molality, and solvent composition.
If no `salts` argument is specified, it will be created via `Clapeyron.auto_binary_salts(model)`
"""
function molality_to_composition(model::ElectrolyteModel,salts::Union{AbstractVector,GroupParam},m,zsolv=SA[1.],ν = salt_stoichiometry(model,salts))
    nc = length(model)
    Z = model.charge
    nions = count(!iszero,Z)
    nneutral = nc - nions
    Mw = mw(model.neutralmodel).*1e-3
    nsalts = salts isa GroupParam ? length(salts.components) : length(salts)
    isalts = 1:nsalts
    iions = 1:nions
    ineutral = 1:nneutral
    if length(zsolv) != nneutral
        throw(error("Incorrect length of zsolv vector,expected length(zsolv) = $nneutral, got $zsolv"))
    end
    ∑mν = sum(m[k]*sum(ν[k,i] for i ∈ iions) for k ∈ isalts)
    ∑zsolv = sum(zsolv)
    TT = Base.promote_eltype(model,m,zsolv)
    x = zeros(TT,length(model))
    ∑xsolvMw = sum(zsolv[j]*Mw[j] for j in ineutral)/∑zsolv
    iion = 0
    for i in 1:length(model)
        if iszero(Z[i])
            x[i] = zsolv[i]/(∑zsolv*(1+∑xsolvMw*∑mν))
        else
            iion += 1
            x[i] = sum(m[k]*ν[k,iion] for k ∈ isalts) / (1/∑xsolvMw+∑mν)
        end
    end
    return x
    #x_solv = zsolv ./ (1+∑xsolvMw*∑mν) ./ ∑zsolv
    #x_ions = [sum(m[k]*ν[k,l] for k ∈ isalts) / (1/∑xsolvMw+∑mν) for l ∈ iions]

    #return vcat(x_solv,x_ions)
end

function molality_to_composition(model::ElectrolyteModel,m,zsolv=SA[1.0],ν = nothing)
    salts = auto_binary_salts(model)
    if isnothing(ν)
        return molality_to_composition(model,salts,m,zsolv)
    else
        return molality_to_composition(model,salts,m,zsolv,ν)
    end
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


auto_binary_salts(model) = auto_binary_salts(model.charge,component_list(model))

function auto_binary_salts(Z,comps)
    #Z = model.charge
    n_ions = count(!iszero,Z)
    res = Tuple{String,Vector{Pair{String,Int}}}[]
    n_ions == 1 && throw(DomainError("cannot create salts with only one ion"))
    n_ions == 0 && return res
    Z_minus = findall(<(0),Z)
    Z_plus = findall(>(0),Z)
    k = 0 #n ions form n-1 independent salts
    for i in 1:length(Z_plus)
        Zi = Z[Z_plus[i]]
        comp_i = comps[Z_plus[i]]
        k == (n_ions - 1) && break
        for j in 1:length(Z_minus)
            k == (n_ions - 1) && break
            k += 1
            Zj = Z[Z_minus[j]]
            comp_j = comps[Z_minus[j]]
            Zij = lcm(abs(Zi),abs(Zj))
            ci = div(Zij,abs(Zi))
            cj = div(Zij,abs(Zj))
            cci = isone(ci) ? "" : string(ci)
            ccj = isone(cj) ? "" : string(cj)
            salt_ij = cci * comp_i * "." *  ccj * comp_j
            push!(res,(salt_ij,[comp_i => ci,comp_j => cj]))
        end
    end
    return res
end

"""
    auto_binary_salts(model)
    auto_binary_salts(Z, components)

Automatically generate a vector of binary salts from ionic components.

# Arguments
- `model`: An electrolyte model containing charge information
- `Z`: Vector of ionic charges for each component
- `components`: Vector of component names

# Returns
A vector of tuples, where each tuple contains:
- `String`: Salt name in the format "n₁Cation.n₂Anion" (e.g., "2Na.SO4" for Na₂SO₄)
- `Vector{Pair{String,Int}}`: Stoichiometric coefficients as component => count pairs

# Description
This function generates `n-1` independent binary salts from `n` ions by pairing cations 
with anions. The stoichiometric coefficients are determined by the least common multiple 
(LCM) of the absolute charges to ensure electroneutrality.

For a system with multiple cations and anions, salts are created by iterating through 
cation-anion pairs until `n-1` salts are formed, which is the number of independent 
salts needed to describe an `n`-ion system.

# Examples
```julia
# For a system with Na⁺ (charge = +1) and Cl⁻ (charge = -1)
Z = [1, -1]
comps = ["Na", "Cl"]
Clapeyron.auto_binary_salts(Z, comps)
# Returns: [("Na.Cl", ["Na" => 1, "Cl" => 1])]

# For a system with Ca²⁺ (charge = +2) and Cl⁻ (charge = -1)
Z = [2, -1]
comps = ["Ca", "Cl"]
Clapeyron.auto_binary_salts(Z, comps)
# Returns: [("Ca.2Cl", ["Ca" => 1, "Cl" => 2])]  # CaCl₂

# For a system with Na⁺, Ca²⁺, and SO₄²⁻ (3 ions → 2 salts)
Z = [1, 2, -2]
comps = ["Na", "Ca", "SO4"]
Clapeyron.auto_binary_salts(Z, comps)
# Returns: [("2Na.SO4", ["Na" => 2, "SO4" => 1]),
#           ("Ca.SO4", ["Ca" => 1, "SO4" => 1])]
```

# Throws
- `DomainError`: If only one ion is present (cannot form salts)

# Notes
- Returns an empty vector if no ions are present (all charges are zero)
- Stoichiometric coefficients of 1 are omitted from salt names (e.g., "Na.Cl" not "1Na.1Cl")
"""
auto_binary_salts

is_electrolyte(model::ElectrolyteModel) = true
 
export molality_to_composition
