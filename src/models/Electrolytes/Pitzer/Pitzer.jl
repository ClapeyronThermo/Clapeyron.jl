abstract type PitzerModel <: ElectrolyteModel end

"""
    PitzerParam <: EoSParam

Parameter container for the Harvie–Møller–Weare / Pitzer ion-interaction
electrolyte activity model (`Pitzer`):
- `charge`: Ionic charge `zᵢ` of each species (`0` for neutral solvents).
- `Mw`: Molar mass `[g·mol⁻¹]` (used to convert solvent mass/molality; `0.01801528 kg·mol⁻¹` for water).
- `beta0`, `beta1`, `beta2`, `Cphi`, `alpha1`, `alpha2`: Cation–anion binary Pitzer interaction parameters.
- `theta`: Like-charge binary mixing parameter `θᵢⱼ` (`cation–cation` or `anion–anion`).
- `psi`: Ternary mixing parameter `ψᵢⱼₖ` for two different like-charged ions `(i, j)` and one opposite-charged ion `k`, stored symmetrically in `(i, j)`.
- `psi_present`: Boolean mask recording which `(i, j, k)` ternary parameters were explicitly supplied in the database or user inputs (even when `0.0`).
- `psi_sources`: Literature source citations parsed from ternary parameter tables.
"""
struct PitzerParam <: EoSParam
    charge::SingleParam{Int}
    Mw::SingleParam{Float64}
    beta0::PairParam{Float64}
    beta1::PairParam{Float64}
    beta2::PairParam{Float64}
    Cphi::PairParam{Float64}
    alpha1::PairParam{Float64}
    alpha2::PairParam{Float64}
    theta::PairParam{Float64}
    psi::Array{Float64,3}
    psi_present::BitArray{3}
    psi_sources::Vector{String}
end

struct Pitzer{I<:IdealModel} <: PitzerModel
    components::Vector{String}
    charge::Vector{Int}
    params::PitzerParam
    idealmodel::I
    temperature_model::Symbol
    require_mixing::Bool
    cm_binary_coeffs::Dict{Tuple{Int,Int,Symbol},NTuple{8,Float64}}
    cm_aphi_coeffs::NTuple{8,Float64}
    references::Vector{String}
end

# Standard Pitzer / Harvie–Møller–Weare constants at 298.15 K
const PITZER_APHI_298 = 0.392
const PITZER_B = 1.2

# Default Christov & Møller (2004) Aϕ(T) correlation coefficients (T in K, 273.15–523.15 K)
const PITZER_DEFAULT_APHI_CM = (
     3.36901532e-01,
    -6.32100430e-04,
     1.92118597e-06,
     0.0,
     9.14252359e+00,
    -1.35143986e-02,
     2.26089488e-03,
     4.52586464e+01
)

default_locations(::Type{Pitzer}) = [
    "properties/molarmass.csv",
    "Electrolytes/properties/charges.csv",
    "Electrolytes/Pitzer/Pitzer_like.csv",
    "Electrolytes/Pitzer/Pitzer_binary.csv",
    "Electrolytes/Pitzer/Pitzer_theta.csv"
]

default_psi_locations(::Type{Pitzer}) = [
    "Electrolytes/Pitzer/Pitzer_psi.csv"
]

default_temperature_locations(::Type{Pitzer}) = [
    "Electrolytes/Pitzer/Pitzer_temperature_ChristovMoller.csv"
]

"""
    Pitzer(solvents, ions;
        idealmodel = BasicIdeal,
        charge = nothing,
        userlocations = String[],
        psi_userlocations = String[],
        temperature_userlocations = String[],
        temperature_model = :constant,
        require_mixing = false,
        verbose = false,
        reference_state = nothing)

## Input parameters
- `charge`: Single Parameter (`Int`) - Ionic charge (`0` for solvent)
- `Mw`: Single Parameter (`Float64`) - Solvent molecular weight `[g·mol⁻¹]`
- `beta0`, `beta1`, `beta2`, `Cphi`, `alpha1`, `alpha2`: Pair Parameters (`Float64`) - Cation–anion binary Pitzer parameters
- `theta`: Pair Parameter (`Float64`) - Like-charge (`cation–cation` or `anion–anion`) mixing parameter `θᵢⱼ`
- `psi`: Ternary Mixing Parameter (`Float64`) - `ψᵢⱼₖ` read from `Pitzer_psi.csv` or `psi_userlocations`

## Description
Harvie–Møller–Weare (HMW) / Pitzer ion-interaction activity model for single
and multicomponent aqueous electrolyte solutions, reusing and extending the
project's validated single-salt (`PitzerSingleSalt`), multi-salt (`PitzerMultiSalt`),
and Christov–Møller temperature-dependent formulations.

Supports:
- Single electrolytes and arbitrary multicomponent salt mixtures (`β⁽⁰⁾, β⁽¹⁾, β⁽²⁾, C^ϕ`, `θᵢⱼ` with Pitzer's unsymmetrical electrostatic mixing `ᴱθᵢⱼ(I)`, and `ψᵢⱼₖ`).
- `temperature_model = :constant` (default `298.15 K` parameterization with `A^ϕ = 0.392`, `b = 1.2`) or `temperature_model = :christov_moller` (`Christov & Møller, 2004` `P(T)` correlations).
- Explicit missing-parameter diagnostics (`pitzer_missing_binary_parameters`, `pitzer_missing_mixing_parameters`) and optional strict enforcement via `require_mixing = true`.

## Examples
```julia
model = Pitzer(["water"], ["sodium", "potassium", "chloride"])
salts = [("NaCl", ["sodium" => 1, "chloride" => 1]),
         ("KCl",  ["potassium" => 1, "chloride" => 1])]
ϕ  = osmotic_coefficient(model, salts, 1e5, 298.15, [1.0, 1.0])
γ± = mean_ionic_activity_coefficient(model, salts, 1e5, 298.15, [1.0, 1.0])
aw = water_activity(model, salts, 1e5, 298.15, [1.0, 1.0])
```

## References
1. Pitzer, K. S. (1973). *Journal of Physical Chemistry*, 77(2), 268–277. [doi:10.1021/j100621a026](https://doi.org/10.1021/j100621a026)
2. Pitzer, K. S., & Kim, J. J. (1974). *Journal of the American Chemical Society*, 96(18), 5701–5707. [doi:10.1021/ja00825a004](https://doi.org/10.1021/ja00825a004)
3. Harvie, C. E., Møller, N., & Weare, J. H. (1984). *Geochimica et Cosmochimica Acta*, 48(4), 723–751. [doi:10.1016/0016-7037(84)90098-X](https://doi.org/10.1016/0016-7037(84)90098-X)
4. Christov, C., & Møller, N. (2004). *Geochimica et Cosmochimica Acta*, 68(6), 1309–1331. [doi:10.1016/j.gca.2003.08.017](https://doi.org/10.1016/j.gca.2003.08.017)
"""
Pitzer

function _match_component_index(cell_str::AbstractString, norm_comps::Vector{String})
    for alias in eachsplit(cell_str, "~|~")
        na = normalisestring(String(alias))
        idx = findfirst(==(na), norm_comps)
        idx !== nothing && return idx
    end
    return nothing
end

function _is_ternary_pitzer_csv(path_or_csv::AbstractString)
    lines = if is_inline_csv(path_or_csv)
        collect(Iterators.take(eachline(IOBuffer(path_or_csv)), 4))
    elseif isfile(path_or_csv)
        open(path_or_csv, "r") do io
            collect(Iterators.take(eachline(io), 4))
        end
    else
        String[]
    end
    length(lines) >= 2 && occursin("ternary", lowercase(lines[2])) && return true
    length(lines) >= 3 && occursin("species3", lowercase(lines[3])) && return true
    return false
end

function _assign_psi_entry!(
    psi::Array{Float64,3},
    psi_present::BitArray{3},
    Z::AbstractVector{<:Integer},
    i1::Int,
    i2::Int,
    i3::Int,
    val::Real
)
    z1, z2, z3 = Z[i1], Z[i2], Z[i3]
    (z1 == 0 || z2 == 0 || z3 == 0) && throw(ArgumentError("Pitzer ternary parameter psi requires three charged ions."))
    like1, like2, opp = if sign(z1) == sign(z2) && sign(z1) != sign(z3)
        (i1, i2, i3)
    elseif sign(z1) == sign(z3) && sign(z1) != sign(z2)
        (i1, i3, i2)
    elseif sign(z2) == sign(z3) && sign(z2) != sign(z1)
        (i2, i3, i1)
    else
        throw(ArgumentError("Pitzer ternary parameter psi requires two like-charged ions and one oppositely charged ion."))
    end
    like1 == like2 && throw(ArgumentError("Pitzer ternary parameter psi requires two distinct like-charged ions."))
    fval = Float64(val)
    psi[like1, like2, opp] = fval
    psi[like2, like1, opp] = fval
    psi_present[like1, like2, opp] = true
    psi_present[like2, like1, opp] = true
    return nothing
end

function _open_pitzer_csv(csv_source::AbstractString)
    if CSV_V1
        return is_inline_csv(csv_source) ?
            CSV.File(IOBuffer(csv_source); header = 3, on_error = :collect, buffer_in_memory = true) :
            CSV.File(csv_source; header = 3, on_error = :collect)
    else
        return is_inline_csv(csv_source) ?
            CSV.File(IOBuffer(csv_source); header = 3, silencewarnings = true, buffer_in_memory = true) :
            CSV.File(csv_source; header = 3, silencewarnings = true)
    end
end

function _read_pitzer_psi_csv!(
    psi::Array{Float64,3},
    psi_present::BitArray{3},
    psi_sources::Vector{String},
    formatted_components::Vector{String},
    Z::AbstractVector{<:Integer},
    csv_source::AbstractString
)
    norm_comps = normalisestring.(formatted_components)
    tbl = _open_pitzer_csv(csv_source)
    cols = Tables.columnnames(tbl)
    norm_col_map = Dict(Symbol(normalisestring(String(c))) => c for c in cols)
    c_sp1 = get(norm_col_map, :species1, nothing)
    c_sp2 = get(norm_col_map, :species2, nothing)
    c_sp3 = get(norm_col_map, :species3, nothing)
    c_psi = get(norm_col_map, :psi, nothing)
    c_src = get(norm_col_map, :source, nothing)
    (c_sp1 === nothing || c_sp2 === nothing || c_sp3 === nothing || c_psi === nothing) && return nothing

    for row in tbl
        s1 = getproperty(row, c_sp1)
        s2 = getproperty(row, c_sp2)
        s3 = getproperty(row, c_sp3)
        pval = getproperty(row, c_psi)
        (ismissing(s1) || ismissing(s2) || ismissing(s3) || ismissing(pval)) && continue
        i1 = _match_component_index(String(s1), norm_comps)
        i2 = _match_component_index(String(s2), norm_comps)
        i3 = _match_component_index(String(s3), norm_comps)
        (i1 === nothing || i2 === nothing || i3 === nothing) && continue
        _assign_psi_entry!(psi, psi_present, Z, i1, i2, i3, pval)
        if c_src !== nothing
            src = getproperty(row, c_src)
            if !ismissing(src) && !isempty(strip(String(src)))
                s_str = String(src)
                s_str in psi_sources || push!(psi_sources, s_str)
            end
        end
    end
    return nothing
end

function _read_pitzer_temperature_csv!(
    cm_binary::Dict{Tuple{Int,Int,Symbol},NTuple{8,Float64}},
    aphi_ref::Base.RefValue{NTuple{8,Float64}},
    formatted_components::Vector{String},
    csv_source::AbstractString
)
    norm_comps = normalisestring.(formatted_components)
    tbl = _open_pitzer_csv(csv_source)
    cols = Tables.columnnames(tbl)
    norm_col_map = Dict(Symbol(normalisestring(String(c))) => c for c in cols)
    c_param = get(norm_col_map, :parameter, nothing)
    c_sp1   = get(norm_col_map, :species1, nothing)
    c_sp2   = get(norm_col_map, :species2, nothing)
    c_a = [get(norm_col_map, Symbol("a$k"), nothing) for k in 1:8]
    (c_param === nothing || any(isnothing, c_a)) && return nothing

    for row in tbl
        pname_raw = getproperty(row, c_param)
        ismissing(pname_raw) && continue
        pname = Symbol(strip(String(pname_raw)))
        coeffs = ntuple(k -> begin
            v = getproperty(row, c_a[k])
            ismissing(v) ? 0.0 : Float64(v)
        end, 8)

        if lowercase(String(pname)) == "aphi"
            aphi_ref[] = coeffs
        elseif c_sp1 !== nothing && c_sp2 !== nothing
            s1 = getproperty(row, c_sp1)
            s2 = getproperty(row, c_sp2)
            (ismissing(s1) || ismissing(s2)) && continue
            i1 = _match_component_index(String(s1), norm_comps)
            i2 = _match_component_index(String(s2), norm_comps)
            (i1 === nothing || i2 === nothing) && continue
            cm_binary[(i1, i2, pname)] = coeffs
            cm_binary[(i2, i1, pname)] = coeffs
        end
    end
    return nothing
end

function Pitzer(
    solvents,
    ions;
    idealmodel = BasicIdeal,
    charge = nothing,
    userlocations = String[],
    psi_userlocations = String[],
    temperature_userlocations = String[],
    temperature_model::Symbol = :constant,
    require_mixing::Bool = false,
    verbose::Bool = false,
    reference_state = nothing
)
    temperature_model in (:constant, :christov_moller) || throw(ArgumentError(
        "Unsupported Pitzer temperature_model = $temperature_model. Expected :constant or :christov_moller."
    ))

    solv_vec = format_components(solvents)
    ion_vec = format_components(ions)
    raw_components = vcat(solv_vec, ion_vec)
    formatted_components = format_components(raw_components)
    nc = length(formatted_components)

    # Separate any ternary psi user inputs from standard 1D/2D getparams userlocations
    extra_psi_files = String[]
    user_psi_dict = nothing
    skip_default_psi = false

    filtered_userlocations = if userlocations isa AbstractVector || userlocations isa Tuple
        std_locs = String[]
        for loc in userlocations
            s = String(loc)
            if s == "@REMOVEDEFAULTS"
                skip_default_psi = true
                push!(std_locs, s)
            elseif is_inline_csv(s) && _is_ternary_pitzer_csv(s)
                push!(extra_psi_files, s)
            elseif !startswith(s, "@") && (isfile(s) || isdir(s))
                for p in getpaths(s)
                    if _is_ternary_pitzer_csv(p)
                        push!(extra_psi_files, p)
                    else
                        push!(std_locs, p)
                    end
                end
            else
                push!(std_locs, s)
            end
        end
        std_locs
    elseif userlocations isa NamedTuple || userlocations isa AbstractDict
        nt_pairs = Pair{Symbol,Any}[]
        for (k, v) in pairs(userlocations)
            ks = Symbol(k)
            if ks == :psi
                user_psi_dict = v
            else
                push!(nt_pairs, ks => v)
            end
        end
        isempty(nt_pairs) ? String[] : (; nt_pairs...)
    else
        userlocations
    end

    Z_vec = Vector{Int}(init_charge(formatted_components, charge; verbose = verbose))
    any(iszero, Z_vec) || throw(ArgumentError("Pitzer requires at least one neutral solvent species (e.g. [\"water\"])."))
    any(>(0), Z_vec)   || throw(ArgumentError("Pitzer requires at least one positively charged cation."))
    any(<(0), Z_vec)   || throw(ArgumentError("Pitzer requires at least one negatively charged anion."))

    params_dict = getparams(
        formatted_components,
        default_locations(Pitzer);
        userlocations = filtered_userlocations,
        ignore_missing_singleparams = ["Mw", "charge", "beta0", "beta1", "beta2", "Cphi", "alpha1", "alpha2", "theta"],
        return_sites = false,
        verbose = verbose
    )

    charge_param = SingleParam("charge", formatted_components, Z_vec)
    Mw_param     = get(params_dict, "Mw", SingleParam("Mw", formatted_components))
    for idx in eachindex(formatted_components)
        if iszero(Z_vec[idx]) && iszero(Mw_param.values[idx])
            throw(ArgumentError("Missing molar mass Mw for solvent '$(formatted_components[idx])'."))
        end
    end

    beta0  = get(params_dict, "beta0",  PairParam("beta0",  formatted_components))
    beta1  = get(params_dict, "beta1",  PairParam("beta1",  formatted_components))
    beta2  = get(params_dict, "beta2",  PairParam("beta2",  formatted_components))
    Cphi   = get(params_dict, "Cphi",   PairParam("Cphi",   formatted_components))
    alpha1 = get(params_dict, "alpha1", PairParam("alpha1", formatted_components))
    alpha2 = get(params_dict, "alpha2", PairParam("alpha2", formatted_components))
    theta  = get(params_dict, "theta",  PairParam("theta",  formatted_components))

    # Populate standard HMW default exponents (alpha1 = 2.0 or 1.4 for 2:2+, alpha2 = 12.0)
    # whenever binary parameters were supplied without explicit alpha1/alpha2 columns.
    for i in 1:nc, j in (i+1):nc
        if (Z_vec[i] > 0 && Z_vec[j] < 0) || (Z_vec[i] < 0 && Z_vec[j] > 0)
            has_bin = !beta0.ismissingvalues[i, j] || !beta1.ismissingvalues[i, j] ||
                      !beta2.ismissingvalues[i, j] || !Cphi.ismissingvalues[i, j]
            if has_bin && alpha1.ismissingvalues[i, j] && iszero(alpha1.values[i, j])
                def_a1 = (abs(Z_vec[i]) >= 2 && abs(Z_vec[j]) >= 2) ? 1.4 : 2.0
                alpha1.values[i, j] = def_a1
                alpha1.values[j, i] = def_a1
                alpha1.ismissingvalues[i, j] = false
                alpha1.ismissingvalues[j, i] = false
            end
            if has_bin && alpha2.ismissingvalues[i, j] && iszero(alpha2.values[i, j])
                alpha2.values[i, j] = 12.0
                alpha2.values[j, i] = 12.0
                alpha2.ismissingvalues[i, j] = false
                alpha2.ismissingvalues[j, i] = false
            end
        end
    end

    psi = zeros(Float64, nc, nc, nc)
    psi_present = falses(nc, nc, nc)
    psi_sources = String[]

    if !skip_default_psi
        for rel_loc in default_psi_locations(Pitzer)
            for p in getpaths(rel_loc; relativetodatabase = true)
                _read_pitzer_psi_csv!(psi, psi_present, psi_sources, formatted_components, Z_vec, p)
            end
        end
    end

    for loc in extra_psi_files
        _read_pitzer_psi_csv!(psi, psi_present, psi_sources, formatted_components, Z_vec, loc)
    end
    for loc in normalize_userlocations(psi_userlocations)
        if is_inline_csv(loc)
            _read_pitzer_psi_csv!(psi, psi_present, psi_sources, formatted_components, Z_vec, loc)
        else
            for p in getpaths(loc)
                _read_pitzer_psi_csv!(psi, psi_present, psi_sources, formatted_components, Z_vec, p)
            end
        end
    end

    if user_psi_dict !== nothing
        if user_psi_dict isa AbstractArray{<:Real,3}
            size(user_psi_dict) == (nc, nc, nc) || throw(DimensionMismatch("psi array must have size ($nc, $nc, $nc)."))
            for i in 1:nc, j in (i+1):nc, k in 1:nc
                v = user_psi_dict[i, j, k]
                if !iszero(v)
                    _assign_psi_entry!(psi, psi_present, Z_vec, i, j, k, v)
                end
            end
        elseif user_psi_dict isa AbstractDict || user_psi_dict isa AbstractVector
            norm_comps = normalisestring.(formatted_components)
            for (k, v) in user_psi_dict
                s1, s2, s3 = k
                i1 = _match_component_index(String(s1), norm_comps)
                i2 = _match_component_index(String(s2), norm_comps)
                i3 = _match_component_index(String(s3), norm_comps)
                (i1 === nothing || i2 === nothing || i3 === nothing) &&
                    throw(ArgumentError("Species in psi key $k not found in model components $formatted_components."))
                _assign_psi_entry!(psi, psi_present, Z_vec, i1, i2, i3, v)
            end
        else
            throw(ArgumentError("Unsupported psi userlocation type $(typeof(user_psi_dict))."))
        end
    end

    cm_binary = Dict{Tuple{Int,Int,Symbol},NTuple{8,Float64}}()
    aphi_ref = Ref{NTuple{8,Float64}}(PITZER_DEFAULT_APHI_CM)
    for rel_loc in default_temperature_locations(Pitzer)
        for p in getpaths(rel_loc; relativetodatabase = true)
            _read_pitzer_temperature_csv!(cm_binary, aphi_ref, formatted_components, p)
        end
    end
    for loc in normalize_userlocations(temperature_userlocations)
        if is_inline_csv(loc)
            _read_pitzer_temperature_csv!(cm_binary, aphi_ref, formatted_components, loc)
        else
            for p in getpaths(loc)
                _read_pitzer_temperature_csv!(cm_binary, aphi_ref, formatted_components, p)
            end
        end
    end

    init_ideal = init_model(idealmodel, raw_components, String[], verbose)
    packaged = PitzerParam(charge_param, Mw_param, beta0, beta1, beta2, Cphi, alpha1, alpha2, theta, psi, psi_present, psi_sources)
    refs = String[
        "10.1021/j100621a026",
        "10.1021/ja00825a004",
        "10.1016/0016-7037(84)90098-X",
        "10.1016/j.gca.2003.08.017"
    ]
    model = Pitzer(formatted_components, Z_vec, packaged, init_ideal, temperature_model, require_mixing, cm_binary, aphi_ref[], refs)
    set_reference_state!(model, reference_state; verbose = verbose)
    return model
end

function Base.show(io::IO, ::MIME"text/plain", model::PitzerModel)
    print(io, "Pitzer Electrolyte Activity Model")
    length(model) == 1 && println(io, " with 1 component:")
    length(model) > 1 && println(io, " with ", length(model), " components:")
    show_comps_with_charge(io, model.components, model.charge)
    print(io, '\n', "Temperature Model: ", model.temperature_model)
    print(io, '\n', "Contains parameters: beta0, beta1, beta2, Cphi, alpha1, alpha2, theta, psi")
    show_reference_state(io, model; space = true)
end

mw(model::PitzerModel) = model.params.Mw.values

"""
    set_pitzer_binary!(model::PitzerModel, cation, anion;
        beta0, beta1, beta2 = 0.0, Cphi = 0.0, alpha1 = 2.0, alpha2 = 12.0)

Register or override binary cation–anion Pitzer parameters in `model`, matching
`set_binary!` from the team's `PitzerMultiSalt` module.
"""
function set_pitzer_binary!(
    model::PitzerModel,
    sp1,
    sp2;
    beta0::Real,
    beta1::Real,
    beta2::Real = 0.0,
    Cphi::Real = 0.0,
    alpha1::Real = 2.0,
    alpha2::Real = 12.0
)
    norm_comps = normalisestring.(model.components)
    i = _match_component_index(String(sp1), norm_comps)
    j = _match_component_index(String(sp2), norm_comps)
    (i === nothing || j === nothing) && throw(ArgumentError("Species ($sp1, $sp2) not found in $(model.components)."))
    sign(model.charge[i]) * sign(model.charge[j]) == -1 ||
        throw(ArgumentError("set_pitzer_binary! requires one cation and one anion."))

    for (par, val) in (
        (model.params.beta0,  Float64(beta0)),
        (model.params.beta1,  Float64(beta1)),
        (model.params.beta2,  Float64(beta2)),
        (model.params.Cphi,   Float64(Cphi)),
        (model.params.alpha1, Float64(alpha1)),
        (model.params.alpha2, Float64(alpha2)),
    )
        par.values[i, j] = val
        par.values[j, i] = val
        par.ismissingvalues[i, j] = false
        par.ismissingvalues[j, i] = false
    end
    return model
end

"""
    set_pitzer_theta!(model::PitzerModel, ion1, ion2, theta_val)

Register or override a like-charge binary mixing parameter `θᵢⱼ` in `model`,
matching `set_theta!` from the team's `PitzerMultiSalt` module.
"""
function set_pitzer_theta!(model::PitzerModel, ion1, ion2, theta_val::Real)
    norm_comps = normalisestring.(model.components)
    i = _match_component_index(String(ion1), norm_comps)
    j = _match_component_index(String(ion2), norm_comps)
    (i === nothing || j === nothing) && throw(ArgumentError("Species ($ion1, $ion2) not found in $(model.components)."))
    i == j && throw(ArgumentError("theta requires two distinct like-charged ions."))
    (model.charge[i] != 0 && sign(model.charge[i]) == sign(model.charge[j])) ||
        throw(ArgumentError("theta requires two ions with the same charge sign."))
    fval = Float64(theta_val)
    model.params.theta.values[i, j] = fval
    model.params.theta.values[j, i] = fval
    model.params.theta.ismissingvalues[i, j] = false
    model.params.theta.ismissingvalues[j, i] = false
    return model
end

"""
    set_pitzer_psi!(model::PitzerModel, ion1, ion2, ion3, psi_val)

Register or override a ternary mixing parameter `ψᵢⱼₖ` in `model`, matching
`set_psi!` from the team's `PitzerMultiSalt` module.
"""
function set_pitzer_psi!(model::PitzerModel, ion1, ion2, ion3, psi_val::Real)
    norm_comps = normalisestring.(model.components)
    i1 = _match_component_index(String(ion1), norm_comps)
    i2 = _match_component_index(String(ion2), norm_comps)
    i3 = _match_component_index(String(ion3), norm_comps)
    (i1 === nothing || i2 === nothing || i3 === nothing) &&
        throw(ArgumentError("Species ($ion1, $ion2, $ion3) not found in $(model.components)."))
    _assign_psi_entry!(model.params.psi, model.params.psi_present, model.charge, i1, i2, i3, psi_val)
    return model
end

# =============================================================================
# CORE VALIDATED PITZER EQUATIONS
# Directly reused and extended from the team's validated `PitzerMultiSalt` and
# `PitzerSingleSalt` (Christov & Møller temperature correlation) modules.
# =============================================================================

"""
    christov_moller_temperature_function(T, a1, a2, a3, a4, a5, a6, a7, a8)

Christov & Møller (2004) 8-term temperature correlation with `T` in Kelvin:
`P(T) = a₁ + a₂T + a₃T² + a₄T³ + a₅/T + a₆ ln(T) + a₇/(T - 263) + a₈/(680 - T)`.
"""
function christov_moller_temperature_function(T, a1, a2, a3, a4, a5, a6, a7, a8)
    return a1 +
           a2 * T +
           a3 * T^2 +
           a4 * T^3 +
           a5 / T +
           a6 * log(T) +
           a7 / (T - 263.0) +
           a8 / (680.0 - T)
end

function _pitzer_Aphi(model::PitzerModel, T)
    if model.temperature_model == :christov_moller
        return christov_moller_temperature_function(T, model.cm_aphi_coeffs...)
    else
        return PITZER_APHI_298 * one(T)
    end
end

function _pitzer_pair_params_at_T(model::PitzerModel, ic::Int, ia::Int, T)
    p = model.params
    a1 = p.alpha1.values[ic, ia]
    a2 = p.alpha2.values[ic, ia]
    if model.temperature_model == :christov_moller && haskey(model.cm_binary_coeffs, (ic, ia, :beta0))
        b0 = christov_moller_temperature_function(T, model.cm_binary_coeffs[(ic, ia, :beta0)]...)
        b1 = haskey(model.cm_binary_coeffs, (ic, ia, :beta1)) ?
            christov_moller_temperature_function(T, model.cm_binary_coeffs[(ic, ia, :beta1)]...) :
            p.beta1.values[ic, ia] * one(T)
        b2 = haskey(model.cm_binary_coeffs, (ic, ia, :beta2)) ?
            christov_moller_temperature_function(T, model.cm_binary_coeffs[(ic, ia, :beta2)]...) :
            p.beta2.values[ic, ia] * one(T)
        cphi = haskey(model.cm_binary_coeffs, (ic, ia, :Cphi)) ?
            christov_moller_temperature_function(T, model.cm_binary_coeffs[(ic, ia, :Cphi)]...) :
            p.Cphi.values[ic, ia] * one(T)
        return b0, b1, b2, cphi, a1, a2
    else
        _1 = one(T)
        return p.beta0.values[ic, ia] * _1,
               p.beta1.values[ic, ia] * _1,
               p.beta2.values[ic, ia] * _1,
               p.Cphi.values[ic, ia]  * _1,
               a1,
               a2
    end
end

function _pitzer_g(x)
    if abs(primalval(x)) < 1e-6
        # Exact Taylor expansion of 2*(1 - (1 + x)*exp(-x))/x^2 around x = 0
        return evalpoly(x, (1.0, -2.0/3.0, 0.25, -1.0/15.0, 1.0/72.0))
    end
    return 2.0 * (1.0 - (1.0 + x) * exp(-x)) / (x * x)
end

# HMW Pitzer g'(x) function: B' is obtained by dividing the beta-weighted sum by I.
function _pitzer_gprime_hmw(x)
    if abs(primalval(x)) < 1e-6
        # Exact Taylor expansion of -2*(1 - (1 + x + 0.5*x^2)*exp(-x))/x^2 around x = 0
        return evalpoly(x, (0.0, -1.0/3.0, 0.25, -0.1, 1.0/36.0))
    end
    return -2.0 * (1.0 - (1.0 + x + 0.5 * x * x) * exp(-x)) / (x * x)
end

function _pitzer_binary_terms(beta0, beta1, beta2, Cphi, alpha1, alpha2, I, zc::Int, za::Int)
    sqrtI = sqrt(I)
    Bphi = beta0 +
           beta1 * exp(-alpha1 * sqrtI) +
           beta2 * exp(-alpha2 * sqrtI)

    B = beta0 +
        beta1 * _pitzer_g(alpha1 * sqrtI) +
        beta2 * _pitzer_g(alpha2 * sqrtI)

    Bprime = if iszero(I)
        zero(B)
    else
        (beta1 * _pitzer_gprime_hmw(alpha1 * sqrtI) +
         beta2 * _pitzer_gprime_hmw(alpha2 * sqrtI)) / I
    end

    C = Cphi / (2.0 * sqrt(Float64(abs(zc * za))))
    return (Bphi = Bphi, B = B, Bprime = Bprime, C = C)
end

# Pitzer (1975) practical approximation to the J(x) integral for unsymmetrical
# electrostatic mixing of like-charged ions of unequal valence.
function _pitzer_J(x)
    primalval(x) <= 1e-14 && return zero(x)
    denominator = 4.0 + 4.581 * x^(-0.7237) * exp(-0.0120 * x^0.528)
    return x / denominator
end

function _pitzer_electrostatic_theta(I, zi::Int, zj::Int, Aphi_val)
    (primalval(I) <= 0.0 || zi == zj) && return zero(I + Aphi_val)
    sqrtI = sqrt(I)
    xij = 6.0 * zi * zj * Aphi_val * sqrtI
    xii = 6.0 * zi^2 * Aphi_val * sqrtI
    xjj = 6.0 * zj^2 * Aphi_val * sqrtI
    return (zi * zj / (4.0 * I)) * (_pitzer_J(xij) - 0.5 * _pitzer_J(xii) - 0.5 * _pitzer_J(xjj))
end

function _pitzer_electrostatic_theta_prime(I, zi::Int, zj::Int, Aphi_val)
    (primalval(I) <= 0.0 || zi == zj) && return zero(I + Aphi_val)
    return Solvers.derivative(ii -> _pitzer_electrostatic_theta(ii, zi, zj, Aphi_val), I)
end

function _pitzer_like_terms(theta_ij, zi::Int, zj::Int, I, Aphi_val)
    Etheta = _pitzer_electrostatic_theta(I, zi, zj, Aphi_val)
    Etheta_prime = _pitzer_electrostatic_theta_prime(I, zi, zj, Aphi_val)
    Phi = theta_ij + Etheta
    Phiprime = Etheta_prime
    Phiphi = Phi + I * Phiprime
    return (Phi = Phi, Phiprime = Phiprime, Phiphi = Phiphi)
end

"""
    pitzer_missing_binary_parameters(model::PitzerModel, m_full = nothing)

Return a vector of `(cation_name, anion_name)` pairs that have no registered
binary Pitzer parameters among the active ions in `model`.
"""
function pitzer_missing_binary_parameters(model::PitzerModel, m_full = nothing)
    Z = model.charge
    p = model.params
    icat = [i for i in eachindex(Z) if Z[i] > 0 && (m_full === nothing || primalval(m_full[i]) > 0.0)]
    iani = [i for i in eachindex(Z) if Z[i] < 0 && (m_full === nothing || primalval(m_full[i]) > 0.0)]
    missing_pairs = Tuple{String,String}[]
    for ic in icat, ia in iani
        has_any = !p.beta0.ismissingvalues[ic, ia] ||
                  !p.beta1.ismissingvalues[ic, ia] ||
                  !p.beta2.ismissingvalues[ic, ia] ||
                  !p.Cphi.ismissingvalues[ic, ia] ||
                  haskey(model.cm_binary_coeffs, (ic, ia, :beta0))
        has_any || push!(missing_pairs, (model.components[ic], model.components[ia]))
    end
    return missing_pairs
end

"""
    pitzer_missing_mixing_parameters(model::PitzerModel, m_full = nothing)

Return a named tuple `(theta = ..., psi = ...)` listing like-charge binary `θ`
and ternary `ψ` parameters absent from `model` for the active ions.
"""
function pitzer_missing_mixing_parameters(model::PitzerModel, m_full = nothing)
    Z = model.charge
    p = model.params
    icat = [i for i in eachindex(Z) if Z[i] > 0 && (m_full === nothing || primalval(m_full[i]) > 0.0)]
    iani = [i for i in eachindex(Z) if Z[i] < 0 && (m_full === nothing || primalval(m_full[i]) > 0.0)]
    missing_theta = Tuple{String,String}[]
    missing_psi = NTuple{3,String}[]

    for a in 1:(length(icat)-1), b in (a+1):length(icat)
        c1, c2 = icat[a], icat[b]
        p.theta.ismissingvalues[c1, c2] && push!(missing_theta, (model.components[c1], model.components[c2]))
        for ia in iani
            p.psi_present[c1, c2, ia] || push!(missing_psi, (model.components[c1], model.components[c2], model.components[ia]))
        end
    end

    for a in 1:(length(iani)-1), b in (a+1):length(iani)
        a1, a2 = iani[a], iani[b]
        p.theta.ismissingvalues[a1, a2] && push!(missing_theta, (model.components[a1], model.components[a2]))
        for ic in icat
            p.psi_present[a1, a2, ic] || push!(missing_psi, (model.components[a1], model.components[a2], model.components[ic]))
        end
    end

    return (theta = missing_theta, psi = missing_psi)
end

"""
    pitzer_molal_properties(model::PitzerModel, T, m_full;
        require_mixing = model.require_mixing,
        check_electroneutrality = true)

Evaluate all molality-scale Pitzer quantities (`I`, `Z`, `phi`, `ln_gamma_molal`,
`gamma_molal`, `ln_aw`, `aw`) at temperature `T` `[K]` for component molality
vector `m_full` (aligned with `model.components`; solvent entries ignored).
"""
function pitzer_molal_properties(
    model::PitzerModel,
    T,
    m_full::AbstractVector;
    require_mixing::Bool = model.require_mixing,
    check_electroneutrality::Bool = true
)
    nc = length(model)
    length(m_full) == nc || throw(DimensionMismatch("Expected molality vector of length $nc, got $(length(m_full))."))

    Z = model.charge
    icat_all = [i for i in 1:nc if Z[i] > 0]
    iani_all = [i for i in 1:nc if Z[i] < 0]
    isolv    = [i for i in 1:nc if Z[i] == 0]

    for i in Iterators.flatten((icat_all, iani_all))
        primalval(m_full[i]) < -1e-15 && throw(ArgumentError("Ion molality must be non-negative (got $(m_full[i]) for $(model.components[i]))."))
    end

    charge_residual = sum(m_full[i] * Z[i] for i in 1:nc)
    if check_electroneutrality
        scale = max(1.0, sum(abs(primalval(m_full[i]) * Z[i]) for i in 1:nc))
        abs(primalval(charge_residual)) <= 1e-10 * scale || throw(ArgumentError(
            "Supplied ion molalities are not electroneutral; charge residual = $(primalval(charge_residual)) mol/kg."
        ))
    end

    missing_bin = pitzer_missing_binary_parameters(model, m_full)
    isempty(missing_bin) || throw(ArgumentError(
        "Missing cation–anion binary Pitzer parameter sets: $missing_bin. " *
        "A mixed solution requires binary parameters for every active cation–anion pair."
    ))

    missing_mix = pitzer_missing_mixing_parameters(model, m_full)
    if require_mixing && (!isempty(missing_mix.theta) || !isempty(missing_mix.psi))
        throw(ArgumentError("Missing mixed-electrolyte Pitzer theta/psi parameters: $missing_mix"))
    end

    TT = Base.promote_eltype(model, T, m_full)
    _0 = zero(TT)
    _1 = one(TT)

    total_m = sum(m_full[i] for i in 1:nc if Z[i] != 0; init = _0)
    I_str   = 0.5 * sum(m_full[i] * (Z[i]^2) for i in 1:nc if Z[i] != 0; init = _0)
    Z_sum   = sum(m_full[i] * abs(Z[i]) for i in 1:nc if Z[i] != 0; init = _0)

    ln_gamma_m = zeros(TT, nc)
    if iszero(total_m)
        gamma_m = ones(TT, nc)
        ln_aw = zeros(TT, length(isolv))
        aw = ones(TT, length(isolv))
        return (
            I = _0,
            Z = _0,
            phi = _1,
            ln_gamma_molal = ln_gamma_m,
            gamma_molal = gamma_m,
            ln_aw = ln_aw,
            aw = aw,
            total_ion_molality = _0,
            charge_residual = charge_residual,
            missing_mixing = missing_mix,
            complete_mixing_model = isempty(missing_mix.theta) && isempty(missing_mix.psi)
        )
    end

    Aphi_val = _pitzer_Aphi(model, T)
    b = PITZER_B

    # Cache binary terms for all cation–anion pairs
    Bphi_mat   = zeros(TT, nc, nc)
    B_mat      = zeros(TT, nc, nc)
    Bprime_mat = zeros(TT, nc, nc)
    C_mat      = zeros(TT, nc, nc)

    for ic in icat_all, ia in iani_all
        # Skip unparameterized inactive pairs (when one of the ions has zero molality)
        if primalval(m_full[ic]) == 0.0 && primalval(m_full[ia]) == 0.0 &&
           model.params.beta0.ismissingvalues[ic, ia] && !haskey(model.cm_binary_coeffs, (ic, ia, :beta0))
            continue
        end
        b0, b1, b2, cphi, a1, a2 = _pitzer_pair_params_at_T(model, ic, ia, T)
        bt = _pitzer_binary_terms(b0, b1, b2, cphi, a1, a2, I_str, Z[ic], Z[ia])
        Bphi_mat[ic, ia]   = bt.Bphi
        B_mat[ic, ia]      = bt.B
        Bprime_mat[ic, ia] = bt.Bprime
        C_mat[ic, ia]      = bt.C
    end

    # Cache like-ion terms for cation–cation and anion–anion pairs
    Phi_mat      = zeros(TT, nc, nc)
    Phiprime_mat = zeros(TT, nc, nc)
    Phiphi_mat   = zeros(TT, nc, nc)

    for a in 1:(length(icat_all)-1), b_idx in (a+1):length(icat_all)
        c1, c2 = icat_all[a], icat_all[b_idx]
        lt = _pitzer_like_terms(model.params.theta.values[c1, c2], Z[c1], Z[c2], I_str, Aphi_val)
        Phi_mat[c1, c2] = Phi_mat[c2, c1] = lt.Phi
        Phiprime_mat[c1, c2] = Phiprime_mat[c2, c1] = lt.Phiprime
        Phiphi_mat[c1, c2] = Phiphi_mat[c2, c1] = lt.Phiphi
    end

    for a in 1:(length(iani_all)-1), b_idx in (a+1):length(iani_all)
        a1, a2 = iani_all[a], iani_all[b_idx]
        lt = _pitzer_like_terms(model.params.theta.values[a1, a2], Z[a1], Z[a2], I_str, Aphi_val)
        Phi_mat[a1, a2] = Phi_mat[a2, a1] = lt.Phi
        Phiprime_mat[a1, a2] = Phiprime_mat[a2, a1] = lt.Phiprime
        Phiphi_mat[a1, a2] = Phiphi_mat[a2, a1] = lt.Phiphi
    end

    sqrtI = sqrt(I_str)
    F = -Aphi_val * (
        sqrtI / (1.0 + b * sqrtI) +
        (2.0 / b) * log(1.0 + b * sqrtI)
    )

    for ic in icat_all, ia in iani_all
        F += m_full[ic] * m_full[ia] * Bprime_mat[ic, ia]
    end
    for a in 1:(length(icat_all)-1), b_idx in (a+1):length(icat_all)
        c1, c2 = icat_all[a], icat_all[b_idx]
        F += m_full[c1] * m_full[c2] * Phiprime_mat[c1, c2]
    end
    for a in 1:(length(iani_all)-1), b_idx in (a+1):length(iani_all)
        a1, a2 = iani_all[a], iani_all[b_idx]
        F += m_full[a1] * m_full[a2] * Phiprime_mat[a1, a2]
    end

    sum_mma_C = _0
    for ic in icat_all, ia in iani_all
        sum_mma_C += m_full[ic] * m_full[ia] * C_mat[ic, ia]
    end

    psi_arr = model.params.psi

    # Cation molality-scale activity coefficients ln(γ_M)
    for M in icat_all
        val = (Z[M]^2) * F
        for ia in iani_all
            val += m_full[ia] * (2.0 * B_mat[M, ia] + Z_sum * C_mat[M, ia])
        end
        for ic in icat_all
            ic == M && continue
            ternary = sum(m_full[ia] * psi_arr[M, ic, ia] for ia in iani_all; init = _0)
            val += m_full[ic] * (2.0 * Phi_mat[M, ic] + ternary)
        end
        for a in 1:(length(iani_all)-1), b_idx in (a+1):length(iani_all)
            a1, a2 = iani_all[a], iani_all[b_idx]
            val += m_full[a1] * m_full[a2] * psi_arr[a1, a2, M]
        end
        val += abs(Z[M]) * sum_mma_C
        ln_gamma_m[M] = val
    end

    # Anion molality-scale activity coefficients ln(γ_X)
    for X in iani_all
        val = (Z[X]^2) * F
        for ic in icat_all
            val += m_full[ic] * (2.0 * B_mat[ic, X] + Z_sum * C_mat[ic, X])
        end
        for ia in iani_all
            ia == X && continue
            ternary = sum(m_full[ic] * psi_arr[X, ia, ic] for ic in icat_all; init = _0)
            val += m_full[ia] * (2.0 * Phi_mat[X, ia] + ternary)
        end
        for a in 1:(length(icat_all)-1), b_idx in (a+1):length(icat_all)
            c1, c2 = icat_all[a], icat_all[b_idx]
            val += m_full[c1] * m_full[c2] * psi_arr[c1, c2, X]
        end
        val += abs(Z[X]) * sum_mma_C
        ln_gamma_m[X] = val
    end

    osmotic_sum = -Aphi_val * (I_str * sqrtI) / (1.0 + b * sqrtI)
    for ic in icat_all, ia in iani_all
        osmotic_sum += m_full[ic] * m_full[ia] * (Bphi_mat[ic, ia] + Z_sum * C_mat[ic, ia])
    end
    for a in 1:(length(icat_all)-1), b_idx in (a+1):length(icat_all)
        c1, c2 = icat_all[a], icat_all[b_idx]
        ternary = sum(m_full[ia] * psi_arr[c1, c2, ia] for ia in iani_all; init = _0)
        osmotic_sum += m_full[c1] * m_full[c2] * (Phiphi_mat[c1, c2] + ternary)
    end
    for a in 1:(length(iani_all)-1), b_idx in (a+1):length(iani_all)
        a1, a2 = iani_all[a], iani_all[b_idx]
        ternary = sum(m_full[ic] * psi_arr[a1, a2, ic] for ic in icat_all; init = _0)
        osmotic_sum += m_full[a1] * m_full[a2] * (Phiphi_mat[a1, a2] + ternary)
    end

    phi = 1.0 + 2.0 * osmotic_sum / total_m

    Mw_vals = model.params.Mw.values
    ln_aw = [-Mw_vals[s] * 1e-3 * phi * total_m for s in isolv]
    aw = exp.(ln_aw)
    gamma_m = exp.(ln_gamma_m)

    return (
        I = I_str,
        Z = Z_sum,
        phi = phi,
        ln_gamma_molal = ln_gamma_m,
        gamma_molal = gamma_m,
        ln_aw = ln_aw,
        aw = aw,
        total_ion_molality = total_m,
        charge_residual = charge_residual,
        missing_mixing = missing_mix,
        complete_mixing_model = isempty(missing_mix.theta) && isempty(missing_mix.psi)
    )
end

# =============================================================================
# CLAPEYRON COMPOSITION, ACTIVITY, AND ELECTROLYTE SOLVER INTERFACES
# =============================================================================

function molality_to_composition(
    model::PitzerModel,
    salts::Union{AbstractVector,GroupParam},
    m,
    zsolv = SA[1.0],
    ν = salt_stoichiometry(model, salts)
)
    nc = length(model)
    Z = component_charges(model)
    nions = count(!iszero, Z)
    nneutral = nc - nions
    Mw = mw(model) .* 1e-3
    nsalts = salts isa GroupParam ? length(salts.components) : length(salts)
    isalts = 1:nsalts
    iions = 1:nions
    ineutral = 1:nneutral
    if length(zsolv) != nneutral
        throw(ArgumentError("Incorrect length of zsolv vector, expected length(zsolv) = $nneutral, got $zsolv"))
    end
    ∑mν = sum(m[k] * sum(ν[k, i] for i in iions) for k in isalts)
    ∑zsolv = sum(zsolv)
    TT = Base.promote_eltype(model, m, zsolv)
    x = zeros(TT, nc)
    isolv_idx = findall(iszero, Z)
    ∑xsolvMw = sum(zsolv[j] * Mw[isolv_idx[j]] for j in ineutral) / ∑zsolv
    iion = 0
    isol = 0
    for i in 1:nc
        if iszero(Z[i])
            isol += 1
            x[i] = zsolv[isol] / (∑zsolv * (1 + ∑xsolvMw * ∑mν))
        else
            iion += 1
            x[i] = sum(m[k] * ν[k, iion] for k in isalts) / (1 / ∑xsolvMw + ∑mν)
        end
    end
    return x
end

"""
    composition_to_molality(model::PitzerModel, z)

Convert a mole-number or mole-fraction vector `z` (aligned with `model.components`)
into a full component molality vector `[mol·kg⁻¹]` relative to total solvent mass.
Solvent entries are set to `0.0`.
"""
function composition_to_molality(model::PitzerModel, z::AbstractVector)
    nc = length(model)
    length(z) == nc || throw(DimensionMismatch("Expected composition vector of length $nc, got $(length(z))."))
    Z = model.charge
    Mw = model.params.Mw.values
    TT = Base.promote_eltype(model, z)
    w_solv_kg = sum(z[i] * Mw[i] * 1e-3 for i in 1:nc if iszero(Z[i]); init = zero(TT))
    primalval(w_solv_kg) > 0.0 || throw(ArgumentError("Solvent mole amount must be positive to define molalities."))
    m_full = zeros(TT, nc)
    for i in 1:nc
        if !iszero(Z[i])
            m_full[i] = z[i] / w_solv_kg
        end
    end
    return m_full
end

function _salts_to_molality_vector(
    model::PitzerModel,
    salts::Union{AbstractVector,GroupParam},
    m,
    ν = salt_stoichiometry(model, salts)
)
    nc = length(model)
    Z = model.charge
    nsalts = size(ν, 1)
    TT = Base.promote_eltype(model, m)
    m_full = zeros(TT, nc)
    iion = 0
    for i in 1:nc
        if !iszero(Z[i])
            iion += 1
            m_full[i] = sum(m[s] * ν[s, iion] for s in 1:nsalts)
        end
    end
    return m_full
end

"""
    ionic_strength(model::PitzerModel, salts, m)
    ionic_strength(model::PitzerModel, z::AbstractVector)

Calculate the molality-scale ionic strength `I = ½ ∑ mᵢ zᵢ²` `[mol·kg⁻¹]`.
"""
function ionic_strength(model::PitzerModel, salts::Union{AbstractVector,GroupParam}, m)
    m_full = _salts_to_molality_vector(model, salts, m)
    return 0.5 * sum(m_full[i] * (model.charge[i]^2) for i in 1:length(model) if !iszero(model.charge[i]))
end

function ionic_strength(model::PitzerModel, z::AbstractVector)
    m_full = composition_to_molality(model, z)
    return 0.5 * sum(m_full[i] * (model.charge[i]^2) for i in 1:length(model) if !iszero(model.charge[i]))
end

reference_chemical_potential_type(::PitzerModel) = :zero

function reference_chemical_potential(model::PitzerModel, p, T, reference = :zero; phase = :unknown, threaded = true, vol0 = nothing)
    _0 = zero(Base.promote_eltype(model, p, T))
    return fill(_0, length(model))
end

function each_split_model(param::PitzerParam, I)
    return PitzerParam(
        each_split_model(param.charge, I),
        each_split_model(param.Mw, I),
        each_split_model(param.beta0, I),
        each_split_model(param.beta1, I),
        each_split_model(param.beta2, I),
        each_split_model(param.Cphi, I),
        each_split_model(param.alpha1, I),
        each_split_model(param.alpha2, I),
        each_split_model(param.theta, I),
        param.psi[I, I, I],
        BitArray(param.psi_present[I, I, I]),
         copy(param.psi_sources)
    )
end

function each_split_model(model::Pitzer, I)
    idx_vec = collect(I)
    inv_map = Dict(old_i => new_i for (new_i, old_i) in enumerate(idx_vec))
    sub_cm = Dict{Tuple{Int,Int,Symbol},NTuple{8,Float64}}()
    for ((i1, i2, pname), coeffs) in model.cm_binary_coeffs
        if haskey(inv_map, i1) && haskey(inv_map, i2)
            sub_cm[(inv_map[i1], inv_map[i2], pname)] = coeffs
        end
    end
    return Pitzer(
        each_split_model(model.components, idx_vec),
        each_split_model(model.charge, idx_vec),
        each_split_model(model.params, idx_vec),
        each_split_model(model.idealmodel, idx_vec),
        model.temperature_model,
        model.require_mixing,
        sub_cm,
        model.cm_aphi_coeffs,
        copy(model.references)
    )
end

"""
    lnγ(model::PitzerModel, p, T, z)

Asymmetric mole-fraction log activity coefficients `ln γₖ⁽ˣ⁾` consistent with
Clapeyron's electrolyte convention (`a_solv = γ_solv⁽ˣ⁾ x_solv` and
`γᵢ⁽ᵐ⁾ = γᵢ⁽ˣ⁾ x_solv`). Satisfies exact Gibbs–Duhem / Euler homogeneity with
[`excess_gibbs_free_energy`](@ref).
"""
function lnγ(model::PitzerModel, p, T, z)
    m_full = composition_to_molality(model, z)
    # Disable electroneutrality check here so component-wise automatic
    # differentiation (∇_z Gᴱ) can perturb individual species freely.
    res = pitzer_molal_properties(model, T, m_full; check_electroneutrality = false)
    Z = model.charge
    nc = length(model)
    ∑z = sum(z)
    z_solv = sum(z[i] for i in 1:nc if iszero(Z[i]))
    x_solv = z_solv / ∑z
    ln_x_solv = log(x_solv)

    TT = Base.promote_eltype(res.phi, ln_x_solv)
    out = zeros(TT, nc)
    isol = 0
    for i in 1:nc
        if iszero(Z[i])
            isol += 1
            x_i = z[i] / ∑z
            # For single-solvent aqueous solution x_i == x_solv and a_w = exp(ln_aw[1])
            out[i] = res.ln_aw[isol] + log(z[i] / z_solv) - log(x_i)
        else
            out[i] = res.ln_gamma_molal[i] - ln_x_solv
        end
    end
    return out
end

function activity_coefficient(
    model::PitzerModel,
    p,
    T,
    z = SA[1.0];
    μ_ref = nothing,
    reference = :zero,
    phase = :unknown,
    threaded = true,
    vol0 = nothing
)
    p̄, T̄, z̄ = ustrip(p, pressure), ustrip(T, temperature), uzstrip(model, z)
    return exp.(lnγ(model, p̄, T̄, z̄))
end

function activity(
    model::PitzerModel,
    p,
    T,
    z = SA[1.0];
    μ_ref = nothing,
    reference = :zero,
    phase = :unknown,
    threaded = true,
    vol0 = nothing
)
    p̄, T̄, z̄ = ustrip(p, pressure), ustrip(T, temperature), uzstrip(model, z)
    ∑z = sum(z̄)
    return exp.(lnγ(model, p̄, T̄, z̄)) .* (z̄ ./ ∑z)
end

function activity_coefficient_impl(model::PitzerModel, p, T, z, μ_ref, reference, phase, threaded, vol0)
    return exp.(lnγ(model, p, T, z))
end

function activity_impl(model::PitzerModel, p, T, z, μ_ref, reference, phase, threaded, vol0)
    ∑z = sum(z)
    return exp.(lnγ(model, p, T, z)) .* (z ./ ∑z)
end

function excess_gibbs_free_energy(model::PitzerModel, p, T, z)
    p̄, T̄, z̄ = ustrip(p, pressure), ustrip(T, temperature), uzstrip(model, z)
    ln_γ_x = lnγ(model, p̄, T̄, z̄)
    return Rgas(model) * T̄ * sum(z̄[i] * ln_γ_x[i] for i in eachindex(z̄))
end

function mean_ionic_activity_coefficient(
    model::PitzerModel,
    salts::Union{AbstractVector,GroupParam},
    p,
    T,
    m,
    zsolv = SA[1.0],
    sat = false
)
    ν = salt_stoichiometry(model, salts)
    m_full = _salts_to_molality_vector(model, salts, m, ν)
    res = pitzer_molal_properties(model, T, m_full)
    iions = findall(!iszero, model.charge)
    ln_γ_ions = @view res.ln_gamma_molal[iions]
    nsalts = size(ν, 1)
    TT = eltype(res.ln_gamma_molal)
    γsm = Vector{TT}(undef, nsalts)
    for s in 1:nsalts
        ν_tot = sum(@view ν[s, :])
        ln_γ_pm = sum(ν[s, k] * ln_γ_ions[k] for k in eachindex(iions)) / ν_tot
        γsm[s] = exp(ln_γ_pm)
    end
    return γsm
end

function mean_ionic_activity_coefficient(model::PitzerModel, p::Number, T, m, zsolv = SA[1.0], sat = false)
    return mean_ionic_activity_coefficient(model, auto_binary_salts(model), p, T, m, zsolv, sat)
end

function mean_ionic_activity_coefficient_sat(
    model::PitzerModel,
    salts::Union{AbstractVector,GroupParam},
    T::Number,
    m,
    zsolv = SA[1.0]
)
    return mean_ionic_activity_coefficient(model, salts, zero(T), T, m, zsolv, true)
end

function mean_ionic_activity_coefficient_sat(model::PitzerModel, T::Number, m, zsolv = SA[1.0])
    return mean_ionic_activity_coefficient_sat(model, auto_binary_salts(model), T, m, zsolv)
end

function osmotic_coefficient(
    model::PitzerModel,
    salts::Union{AbstractVector,GroupParam},
    p,
    T,
    m,
    zsolv = SA[1.0],
    sat = false
)
    ν = salt_stoichiometry(model, salts)
    m_full = _salts_to_molality_vector(model, salts, m, ν)
    res = pitzer_molal_properties(model, T, m_full)
    nsolv = count(iszero, model.charge)
    return fill(res.phi, nsolv)
end

function osmotic_coefficient(model::PitzerModel, p, T, m, zsolv = SA[1.0], sat = false)
    return osmotic_coefficient(model, auto_binary_salts(model), p, T, m, zsolv, sat)
end

function osmotic_coefficient_sat(
    model::PitzerModel,
    salts::Union{AbstractVector,GroupParam},
    T,
    m,
    zsolv = SA[1.0]
)
    return osmotic_coefficient(model, salts, zero(T), T, m, zsolv, true)
end

"""
    water_activity(model::PitzerModel, salts, p, T, m, zsolv = SA[1.0])
    water_activity(model::PitzerModel, p, T, m, zsolv = SA[1.0])

Calculate the thermodynamic solvent (water) activity `a_w = exp(-M_w ϕ ∑ mᵢ)`
for `model` at temperature `T` `[K]` and salt molalities `m` `[mol·kg⁻¹]`.
"""
function water_activity(
    model::PitzerModel,
    salts::Union{AbstractVector,GroupParam},
    p,
    T,
    m,
    zsolv = SA[1.0]
)
    ν = salt_stoichiometry(model, salts)
    m_full = _salts_to_molality_vector(model, salts, m, ν)
    res = pitzer_molal_properties(model, T, m_full)
    return length(res.aw) == 1 ? res.aw[1] : res.aw
end

function water_activity(model::PitzerModel, p::Number, T, m, zsolv = SA[1.0])
    return water_activity(model, auto_binary_salts(model), p, T, m, zsolv)
end

"""
    molal_activity_coefficient(model::PitzerModel, salts, p, T, m)
    molal_activity_coefficient(model::PitzerModel, p, T, m)

Return the vector of individual-ion molality-scale activity coefficients `γᵢ⁽ᵐ⁾`
for the charged species in `model`.
"""
function molal_activity_coefficient(
    model::PitzerModel,
    salts::Union{AbstractVector,GroupParam},
    p,
    T,
    m
)
    ν = salt_stoichiometry(model, salts)
    m_full = _salts_to_molality_vector(model, salts, m, ν)
    res = pitzer_molal_properties(model, T, m_full)
    iions = findall(!iszero, model.charge)
    return res.gamma_molal[iions]
end

function molal_activity_coefficient(model::PitzerModel, p::Number, T, m)
    return molal_activity_coefficient(model, auto_binary_salts(model), p, T, m)
end

export Pitzer, PitzerModel,
       set_pitzer_binary!, set_pitzer_theta!, set_pitzer_psi!,
       composition_to_molality, ionic_strength,
       pitzer_missing_binary_parameters, pitzer_missing_mixing_parameters,
       pitzer_molal_properties, water_activity, molal_activity_coefficient,
       christov_moller_temperature_function
