#=
LS: liquid-state theory (Zhang et al. 2016) for a group-contribution
pseudo-component built from charged and neutral beads (e.g. a coarse-
grained polyampholyte/IDP chain, or a free ion species). Every group must
be named "+" (charge +1), "-" (charge -1), or "0" (neutral) -- charge is
read directly off the group name, not supplied separately.

Constructed directly from already-parsed group/bond data, exactly like
`HeterogcPCPSAFT`'s own `input` format, since a synthetic (non-database)
pseudo-component has no CSV entry to look up from:

    LS([("polymer", ["+"=>10, "-"=>10], [("+","-")=>1, ("+","+")=>4, ("-","-")=>4]),
        ("counterion", ["-"=>1], [])])

Each tuple is `(component_name, group=>count pairs, (group1,group2)=>bond_count pairs)`.
Deriving this input from a linear bead *sequence* (e.g. `[1,-1,1,-1,...]`)
is not this model's concern -- that belongs in the caller, the same way
resolving a component *name* into groups (SMILES lookup, a fixed CSV table,
...) is HeterogcPCPSAFT's caller's concern, not HeterogcPCPSAFT itself.

`LS` is a composite of `LSNeutral` (hard-sphere + TPT1 chain term) and
`LSIon` (restricted-primitive-model Blum-MSA electrostatics + the
charge-induced chain-term correction), following the same neutral-model/
ion-model composition pattern as `ESElectrolyte`. It subtypes the bare
`ElectrolyteModel` (not `ESElectrolyteModel`): `LS` supplies its own
`lb_volume`/`T_scale`, and none of `ESElectrolyteModel`'s generic
`mw`/`p_scale`/`x0_volume_liquid`/`x0_volume_gas` fallbacks (which assume a
real-named-species solvent) apply or are needed here.

Reference: Zhang, Shakhnovich, et al. (Zhang, Y.; Shen, K.; Liu, K.; Cook,
S. R.; Wingreen, N. S.; ... 2016).
=#

const LS_SIGMA = 3e-10 # m -- bead diameter, shared by every bead in every component

const _LS_GROUP_CHARGE = Dict("+" => 1, "-" => -1, "0" => 0)

"""
    _ls_group_input(input)

Split the user-facing `(name, groups, intragroups)` triples into the
`(name, groups)` pairs `GroupParam` itself takes, plus a `gc_intragroups`
list ready for `build_gc_intragroups!`. A component with no intragroup
bonds given and exactly one distinct group is filled in with the trivial
zero-bond self-pair (the same convention `GroupParam`'s own CSV-loading
path uses for a single-group component) -- e.g. a bare free ion needs no
bonds at all, but `build_gc_intragroups!` treats a genuinely empty
per-component list as missing data, not "zero bonds".
"""
function _ls_group_input(input)
    group_input = [(name, grouppairs) for (name, grouppairs, _) in input]
    gc_intragroups = Vector{Vector{Pair{Tuple{String,String},Int}}}(undef, length(input))
    for (c, (name, grouppairs, bonds)) in enumerate(input)
        if isempty(bonds)
            gnames = first.(grouppairs)
            length(gnames) == 1 ||
                throw(ArgumentError("component \"$name\" has $(length(gnames)) distinct groups but no intragroup bonds were given"))
            gc_intragroups[c] = [(gnames[1], gnames[1]) => 0]
        else
            gc_intragroups[c] = bonds
        end
    end
    return group_input, gc_intragroups
end

"""
    _ls_group_charge(g::String)

Charge for a group named `g`. 
Reads only the *first character* against the fixed `"+"`/`"-"`/`"0"` convention (see this file's header), not the whole
name, so a caller needing globally-unique group names for other reasons
(e.g. one group per individual bead *position*, `"+_1"`, `"-_2"`, ...) can
suffix freely -- the leading charge symbol is all that's load-bearing.
"""
function _ls_group_charge(g::String)
    isempty(g) && throw(ArgumentError("empty LS group name"))
    symbol = string(first(g))
    haskey(_LS_GROUP_CHARGE, symbol) ||
        throw(ArgumentError("unrecognized LS group name \"$g\" -- must start with \"+\", \"-\", or \"0\""))
    return _LS_GROUP_CHARGE[symbol]
end

function _ls_build_groups(input)
    group_input, gc_intragroups = _ls_group_input(input)
    groups = GroupParam(group_input)
    build_gc_intragroups!(groups, gc_intragroups)

    Zvalues = _ls_group_charge.(groups.flattenedgroups)
    Z = SingleParam("charge", groups.flattenedgroups, Zvalues)
    return groups, Z
end

"""
    _ls_bond_counts(groups, Zvalues, c)

Given a `GroupParam`, a vector of charges per bead and and index `c`, returns:
- `Npm`: number of bonds between unlike-charge groups
- `Nnc` number of bonds involving a neutral group (no units)

The `GroupParam` must have their intragroup matrix initialized and, and only the group names (`+`,`-`,`0`) are allowed.
"""
function _ls_bond_counts(groups::GroupParam, Zvalues::Vector{Int}, c::Int)

    #=
    `groups.n_intergroups[c]` (already sized to the *global* flattened-group
    count by `build_gc_intragroups!`) and the per-group charges `Zvalues`. A
    same-charge bond (`+`-`+` or `-`-`-`, whether between two distinct groups
    or a group's own self-count) contributes to neither; a `0`-`0` self-count
    does contribute to `Nnc`, so the diagonal is handled explicitly alongside
    the off-diagonal pairs.
    =#
    n_mat = groups.n_intergroups[c]
    i_groups_c = groups.i_groups[c]
    Npm = 0
    Nnc = 0
    for idx1 in eachindex(i_groups_c)
        k = i_groups_c[idx1]
        Zvalues[k] == 0 && (Nnc += n_mat[k, k])
        for idx2 in idx1+1:length(i_groups_c)
            l = i_groups_c[idx2]
            n_mat[k, l] == 0 && continue
            if Zvalues[k] == 0 || Zvalues[l] == 0
                Nnc += n_mat[k, l]
            elseif Zvalues[k] != Zvalues[l]
                Npm += n_mat[k, l]
            end
        end
    end
    return Npm, Nnc
end

#=
Per-component `N`/`Npm`/`Nnc`, each a `SingleParam{Int}` indexed by
`components`: `N` is the total bead count, `Npm`/`Nnc` the unlike-charge/
neutral-involving bond counts (see [`_ls_bond_counts`](@ref)). A free ion
component (`N=1,Npm=0,Nnc=0`) makes both `LSNeutral`/`LSIon`'s per-component
chain-term coefficients `(1-N)` and `(1+2Npm+Nnc-N)` vanish identically, so
no special-casing is needed elsewhere for its presence.
=#
struct LSNeutralParam <: EoSParam
    N::SingleParam{Int}   #
    Npm::SingleParam{Int} #`Npm` (bonds between unlike-charge groups)
    Nnc::SingleParam{Int} #`Nnc` (bonds involving a neutral group)
end

struct LSIonParam <: EoSParam
    Z::SingleParam{Int}    # charge per GROUP, aligned with groups.flattenedgroups
    N::SingleParam{Int}
    Npm::SingleParam{Int}
    Nnc::SingleParam{Int}
end

"""
    ls_group_densities(groups, V, z)

Reduced per-group number densities `ρ★[k] = N_A Σ_c z[c] n_flattenedgroups[c][k] σ³ / V`
(`σ`=`LS_SIGMA`), summed over every component `c`.
"""
function ls_group_densities(groups::GroupParam, V, z)
    TT = Base.promote_eltype(groups,V,z)
    ρ = similar(z,TT,length(groups.n_flattenedgroups[1]))
    η = N_A * LS_SIGMA^3 / V
    for c in 1:length(z)
        zc = z[c]
        nc = groups.n_flattenedgroups[c]
        for k in 1:length(ρ)
            ρ[k] = zc * nc[k] * η
        end
    end
    return ρ
end

struct LSNeutral <: EoSModel
    components::Vector{String}
    groups::GroupParam
    params::LSNeutralParam
end

"""
    LSNeutral(input)

## Model Parameters
- `N`: Bead count (no units)
- `Npm`: Bonds between unlike-charge groups (no units)
- `Nnc` Bonds involving a neutral group (no units)

## Input models
- `idealmodel`: Ideal Model

## Description
Neutral (density-only) half of the LS-theory bulk EOS. 
Hard-sphere (BMCSL) + the Γ_MSA=0 limit of the TPT1 chain term. 
See [`LSIon`](@ref) for the charge-induced half; [`LS`](@ref) composes the two and is the usual entry point -- see its docstring for `input`'s format.

`a_res` depends only on `N`/`Npm`/`Nnc` and the total reduced density (invariant to how beads are grouped).
"""
function LSNeutral(input)
    groups, Z = _ls_build_groups(input)
    components = groups.components
    Zvalues = Z.values
    ncomp = length(components)
    Nv, Npmv, Nncv = zeros(Int, ncomp), zeros(Int, ncomp), zeros(Int, ncomp)
    for c in 1:ncomp
        Nv[c] = sum(groups.n_groups[c])
        Npmv[c], Nncv[c] = _ls_bond_counts(groups, Zvalues, c)
    end
    params = LSNeutralParam(SingleParam("N", components, Nv), SingleParam("Npm", components, Npmv), SingleParam("Nnc", components, Nncv))
    return LSNeutral(components, groups, params)
end

function a_res(model::LSNeutral, V, T, z, ρ★ = ls_group_densities(model.groups, V, z))
    p = model.params
    ρtot = sum(ρ★)
    η = (π/6) * ρtot
    fhs = 6 * η^2 * (4 - 3η) / (π * (1 - η)^2)

    # Hard-sphere contact-value correlation matching the classical-DFT chain
    # term's own kernel (the BMCSL-style form, not Zhang et al. 2016's own
    # yhs=(2+η)/(2(1-η)²) closed form -- those differ by a few percent at
    # typical η, so using this one keeps the bulk model and a reused DFT
    # kernel mutually consistent).
    c1 = 1 / (1 - η)
    c2 = 3η / (1 - η)^2
    c3 = 2η^2 / (1 - η)^3
    r = 0.5
    g_hs = c1 + r * c2 + r^2 * c3

    # Γ_MSA=0 limit: every bond's contact value collapses to g_hs, giving
    # (1-N[c])*log(g_hs) per component -- a free (unbonded) ion component
    # (N=1) contributes exactly 0, no special-casing needed.
    ρchain_tot = zero(ρtot)
    fch0 = zero(ρtot)
    Nvals = p.N.values
    for c in eachindex(z)
        ρchain_c = N_A * z[c] * LS_SIGMA^3 / V
        ρchain_tot += ρchain_c
        fch0 += ρchain_c * (1 - Nvals[c]) * log(g_hs)
    end

    # a_res is per mole of the whole mixture (Σ_c z[c] chains/ions), not per
    # mole of monomer/bead of any single component -- fhs/fch0 are energy
    # densities, so the correct normalization is ρchain_tot, not the
    # aggregate bead density ρtot.
    return (fhs + fch0) / ρchain_tot
end

struct LSIon{ϵ} <: EoSModel
    components::Vector{String}
    groups::GroupParam
    params::LSIonParam
    RSPmodel::ϵ
end

"""
    LSIon(input; RSPmodel=ConstRSP())

## Model Parameters
- `Z`: Charge of each bead (no units)
- `N`: Bead count (no units)
- `Npm`: Bonds between unlike-charge groups (no units)
- `Nnc` Bonds involving a neutral group (no units)

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
Charge-induced half of the LS-theory bulk EOS.
The restricted-primitive-model (equal bead size) Blum-MSA electrostatic free energy, plus the charge-induced correction to the TPT1 chain term.
Shares its bead/bond data with [`LSNeutral`](@ref); [`LS`](@ref) composes the two and is the usualn entry point -- see its docstring for `input`'s format.
"""
function LSIon(input; RSPmodel=ConstRSP())
    groups, Z = _ls_build_groups(input)
    components = groups.components
    Zvalues = Z.values
    ncomp = length(components)
    Nv, Npmv, Nncv = zeros(Int, ncomp), zeros(Int, ncomp), zeros(Int, ncomp)
    for c in 1:ncomp
        Nv[c] = sum(groups.n_groups[c])
        Npmv[c], Nncv[c] = _ls_bond_counts(groups, Zvalues, c)
    end
    params = LSIonParam(Z, SingleParam("N", components, Nv), SingleParam("Npm", components, Npmv), SingleParam("Nnc", components, Nncv))
    return LSIon(components, groups, params, RSPmodel)
end

data(model::LSIon, V, T, z) = ls_group_densities(model.groups, V, z)

function a_res(model::LSIon, V, T, z, ρ★ = @f(data))
    Z = model.params.Z.values
    N = model.params.N.values
    Npm = model.params.Npm.values
    Nnc = model.params.Npm.values

    ϵr = dielectric_constant(model.RSPmodel, 1.0, T, z)
    lB = e_c^2 / (4π * ϵ_0 * ϵr * LS_SIGMA * k_B * T)

    # Restricted primitive model (all beads the same size): the Blum-MSA
    # screening parameter Γ_MSA has this closed-form solution -- no
    # fixpoint iteration needed (unlike the general/asymmetric-size MSA.jl).
    ρZ3 = @sum(ρ★[i] * Z[i]*Z[i]*abs(Z[i]))
    κ_MSA = sqrt(4π * lB * ρZ3)
    Γ_MSA = (-1 + sqrt(1 + 2κ_MSA)) / 2
    fel = -Γ_MSA^3 * (2/3 + Γ_MSA) / π

    # log(ypp) = log(yhs) - lB/(1+Γ_MSA)² + lB, log(ypm) = log(yhs) +
    # lB/(1+Γ_MSA)² - lB; substituting into fch(Γ_MSA) and subtracting
    # fch(0), every log(yhs) term cancels exactly, leaving this closed
    # form. Summed per component -- a free ion component contributes
    # exactly 0 here too (its (1+2·0+0-1)=0 identically).
    ρchain_tot = zero(κ_MSA)
    Δfch = zero(κ_MSA)
    
    for c in eachindex(z)
        ρchain_c = N_A * z[c] * LS_SIGMA^3 / V
        ρchain_tot += ρchain_c
        Δfch += ρchain_c * (1 + 2*Npm[c] + Nnc[c] - N[c]) * lB * (1 - 1 / (1 + Γ_MSA)^2)
    end

    return (fel + Δfch) / ρchain_tot
end


#=

. 
`LS` subtypes `Clapeyron.ElectrolyteModel` -- not `ESElectrolyteModel` -- so
none of that type's generic `mw`/`p_scale`/`T_scale`/`lb_volume`/ `x0_volume_liquid`/`x0_volume_gas` fallbacks 
(tuned for a real-named-solvent + dissolved-salt representation) apply; 
`LS` defines its own `lb_volume` and `T_scale` below.

=#
struct LS{N<:EoSModel,I<:EoSModel} <: ElectrolyteModel
    components::Vector{String}
    neutralmodel::N
    ionmodel::I
    idealmodel::BasicIdeal
    groups::GroupParam
    charge::Vector{Int}
end

"""
    LS(input; RSPmodel=ConstRSP())

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
Bulk liquid-state-theory EOS (Zhang et al. 2016) for a group-contribution pseudo-component built from charged/neutral beads. 
`input` is a list of `(component_name, group=>count pairs, (group1,group2)=>bond_count pairs)` triples, exactly the `HeterogcPCPSAFT`-style convention:

```julia
model = LS([("polymer", ["+"=>10, "-"=>10], [("+","-")=>1, ("+","+")=>4, ("-","-")=>4]),
            ("counterion", ["-"=>1], [])])
```

Every group must be named `"+"` (charge +1), `"-"` (charge -1), or `"0"` (neutral) -- charge is read directly off the group name. 
A component with a single distinct group and no bonds given (a bare free ion) has its trivial zero-bond self-pair filled in automatically.

## References
1. Zhang, P., Alsaifi, N. M., Wu, J., & Wang, Z.-G. (2016). Salting-out and salting-in of polyelectrolyte solutions: A liquid-state theory study. Macromolecules, 49(24), 9720–9730. [doi:10.1021/acs.macromol.6b02160](https://doi.org/10.1021/acs.macromol.6b02160)
"""
function LS(input; RSPmodel=ConstRSP())
    neutralmodel = LSNeutral(input)
    ionmodel = LSIon(input; RSPmodel=RSPmodel)
    return LS(neutralmodel.components, neutralmodel, ionmodel, BasicIdeal(), neutralmodel.groups, ionmodel.params.Z.values)
end

data(model::LS, V, T, z) = ls_group_densities(model.groups, V, z)

function a_res(model::LS, V, T, z, ρ★ = @f(data))
    return a_res(model.neutralmodel, V, T, z, ρ★) + a_res(model.ionmodel, V, T, z, ρ★)
end

function lb_volume(model::LS, z)
    Nvals = model.neutralmodel.params.N.values
    return (π/6) * N_A * LS_SIGMA^3 * sum(z[c] * Nvals[c] for c in eachindex(z))
end

function T_scale(model::LS, z)
    ϵr = dielectric_constant(model.ionmodel.RSPmodel, 1.0, 298.15, z)
    # nominal scale: the temperature at which the Bjerrum length equals one
    # bead diameter (lB=1 in the reduced units used throughout this model).
    return e_c^2 / (4π * ϵ_0 * ϵr * LS_SIGMA * k_B)
end

"""
    lB_to_T(model::LS, lB; z=SA[1.0])

Convert a target (reduced) Bjerrum length to the corresponding real
temperature for `model`'s ion submodel/dielectric constant.
"""
function lB_to_T(model::LS, lB::Real; z=SA[1.0])
    ϵr = dielectric_constant(model.ionmodel.RSPmodel, 1.0, 298.15, z)
    return e_c^2 / (4π * ϵ_0 * ϵr * LS_SIGMA * k_B * lB)
end

#=
    x0_crit_pure(model::LS, z=SA[1.0])

LS's critical point is electrostatically driven with a much lower critical packing fraction (~0.005-0.02) than the generic `x0_crit_pure` default (~0.3) assumes, 
which converges to a spurious point here -- delegate to the general bisection-based warm start instead (see `x0_crit_pure_bisection`).
=#
x0_crit_pure(model::LS, z=SA[1.0]) = x0_crit_pure_bisection(model, z)

#=

The generic `x0_crit_mix` default calls `split_pure_model(model)` then `crit_pure` on each isolated pure component -- for LS's 1-group/zero-bond free-ion component, 
no `is_splittable`/`default_splitter` method is defined, so that throws a `BoundsError`. 
Delegate to the general fixed-ray bisection instead (see `x0_crit_mix_bisection`).
=#
x0_crit_mix(model::LS, z) = x0_crit_mix_bisection(model, z)


function component_charges(model::LS)
    return [dot(model.charge,model.groups.n_flattenedgroups[c]) for c in 1:length(model)]
end

export LS, LSNeutral, LSIon
export component_charges
