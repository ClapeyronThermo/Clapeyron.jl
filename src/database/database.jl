@enum CSVType invaliddata namedtupledata singledata pairdata assocdata groupdata
const NO_KIJ = """@REPLACE Clapeyron Database File
no Parameters [csvtype = unlike]
species1,species2,k
"""

const NO_ASSOC = """@REPLACE Clapeyron Database File
no Parameters [csvtype = assoc]
species1,species2,site1,site2,epsilon_assoc,bondvol
"""

include("database_rawparam.jl")
include("database_utils.jl")

"""
    params = getparams(components,locations;kwargs...)
Returns a `Dict{String,ClapeyronParam}` containing all the parameters found for the list of components
in the available CSVs. `locations` are the locations relative to `Clapeyron` database. The available keywords are the ones used ∈ [`ParamOptions`](@ref)
if `return_sites` is set to true, `getparams` will add a "sites" value in the params result, containing a `SiteParam` built with the input parameters.

## Single to Pair promotion

When reading multiple CSVs, if a parameter's name appears in a single parameter file and in a pair parameter file, the single parameter values will be promoted to be the diagonal values of the pair interaction matrix:

**`my_parameter_single.csv`**
```
Clapeyron Database File
like parameters
species,a
sp1,1000
sp2,700
sp3,850
```
**`my_parameter_pair.csv`**
```
Clapeyron Database File
pair parameters
species1,species2,a
sp1,sp2,875
sp2,sp3,792
sp3,sp1,960

julia> res = getparams(["sp1","sp2"],userlocations = [my_parameter_single.csv,my_parameter_pair.csv])
Dict{String, Clapeyron.ClapeyronParam} with 1 entry:
  "a" => PairParam{Int64}("a")["sp1", "sp2"]

julia> res["a"].values
2×2 Matrix{Int64}:
 1000  875
  875  700
```
This promotion is only supported for Single-Pair combinations. Other CSV type combinations will fail.

## In-memory CSV parsing

If you pass any string starting with `Clapeyron Database File`, it will be parsed as a CSV instead of being used as a filepath:

```julia-repl
julia> x = \"\"\"Clapeyron Database File,
       in memory like parameters
       species,a,b
       sp1,1000,0.05
       sp2,700,0.41
       \"\"\"
"Clapeyron Database File,\nin memory parameters [csvtype = like,grouptype = in_memory_read]\nspecies,a,b\nsp1,1000,0.05\nsp2,700,0.41\n"
julia> Clapeyron.getparams(["sp1","sp2"],userlocations = [x])
Dict{String, Clapeyron.ClapeyronParam} with 2 entries:
  "b" => SingleParam{Float64}("b")["sp1", "sp2"]
  "a" => SingleParam{Int64}("a")["sp1", "sp2"]
```
## Special prefixes

There are some special prefixes that are used by the parser to signal some specific behaviour to be done at parsing time, for one CSV or a group of them:
- `@DB`: replaces the path by the current Clapeyron default database. When doing `getparams(components,["location"])`, the paths are lowered to `getparams(components,userlocations = ["@DB/location"])`.
In a way, is a path shortcut used internally by Clapeyron to parse it's own database. You can change the path where `@DB` points to (or add other path shortcuts), via adding a corresponding entry to the `Clapeyron.SHORT_PATHS` Dict.
- `@REPLACE`: Any filepath starting with `@REPLACE` will clear all previous appearances of the parameter names found in the CSV that contains the prefix.
- `@REMOVEDEFAULTS`: it is used alone, and needs to be passed at the first position of the vector of `userlocations`. It will skip parsing of the default parameters:

The effect of the parser can be summarized by the following examples:

```
model = PCSAFT(["water"],userlocations = ["@REMOVEDEFAULTS"]) #fails, no parameters found, no CSV parsed
model = PCSAFT(["water"],userlocations = ["@REPLACE/empty_params.csv"]) #fails, no parameters found, default parameters parsed and then removed
model = PCSAFT(["water"],userlocations = ["@REPLACE/my_pcsaft_kij.csv"]) #success, default kij parameters replaced by the ones on `my_pcsaft_kij.csv`
model = PCSAFT(["water"],userlocations = ["@REMOVEDEFAULTS","@DB/SAFT/PCSAFT","@DB/properties/molarmass.csv"]) #sucess. Default parameters csv removed, and parsed again, using the @DB prefix to point to the default database.
```

You can use the `@REPLACE` keyword in an in-memory CSV by adding it at the start of the string, followed by a space:
```
#This will replace all previous parsed occurences of `a` and `b`
x_replace = \"\"\"@REPLACE Clapeyron Database File,
in memory like parameters
species,a,b
sp1,1000,0.05
sp2,700,0.41
\"\"\"
```

## CSV type detection and group type

The second line of the csv is used for comments and to identify the type of CSV used. For example:

    ```
x = \"\"\"Clapeyron Database File
       in memory like parameters
       species,a,b
       sp1,1000,0.05
       sp2,700,0.41
       \"\"\"
```
Will be parsed as a table with single parameter data. If you want more flexibility, you can instead pass the csvtype between brackets:

```
x = \"\"\"Clapeyron Database File
       i can write anything here, unlike, association [csvtype = like] but the csv type is already specified.
       species,a,b
       sp1,1000,0.05
       sp2,700,0.41
       \"\"\"
```
Additionaly, there are some cases when you want to be absolutely sure that your types don't clash with the default values. This is the case with different group parametrizations of UNIFAC (Dortmund, VTPR, PSRK):

```
julia> model = UNIFAC(["methanol","ethanol"])
UNIFAC{PR{BasicIdeal, PRAlpha, NoTranslation, vdW1fRule}} with 2 components:
 "methanol": "CH3OH" => 1
 "ethanol": "CH2" => 1, "CH3" => 1, "OH (P)" => 1
Group Type: UNIFACDortmund
Contains parameters: A, B, C, R, Q

julia> model = PSRKUNIFAC(["methanol","ethanol"])
UNIFAC{BasicIdeal} with 2 components:
 "methanol": "CH3OH" => 1
 "ethanol": "CH2" => 1, "CH3" => 1, "OH" => 1
Group Type: PSRK
Contains parameters: A, B, C, R, Q
```

The models are the same (`UNIFAC`), but the group parametrizations are different. This is specified with the `grouptype` keyword. For example, if we see `UNIFAC_groups.csv`, it starts with:

```
Clapeyron Database File,
modified UNIFAC (Dortmund) Groups [csvtype = groups,grouptype = UNIFACDortmund]
species,groups
ethane,"[""CH3"" => 2]"
propane,"[""CH3"" => 2, ""CH2"" => 1]"
butane,"[""CH3"" => 2, ""CH2"" => 2]"
...
```

For compatibility reasons, if you pass a CSV without grouptype, it will be accepted, but two CSV with different specified group types cannot be merged:

```
x1 = \"\"\"Clapeyron Database File
       parameterization 1 [csvtype = like,grouptype = param1]
       species,a,b
       sp1,1000,0.05
       sp2,700,0.41
       \"\"\"
x2 = \"\"\"Clapeyron Database File
       fitted to data [csvtype = like,grouptype = fitted]
       species,a,b
       sp1,912,0.067
       sp2,616,0.432
       \"\"\"
```

If we pass the same parameters, with different group types, the parser will fail.

```julia-repl
julia> Clapeyron.getparams(["sp1","sp2"],userlocations = [x1,x2])
ERROR: cannot join two databases with different group types:
current group type: param1
incoming group type: fitted
```

Note, that the parser will not fail if you pass different parameters with different group types (For example if `a` has `param1` group type and `b` has `fit` group type)
"""
function getparams(components,
                    locations=String[];
                    userlocations = String[],
                    asymmetricparams=String[],
                    ignore_missing_singleparams=String[],
                    ignore_headers = IGNORE_HEADERS,
                    verbose::Bool=false,
                    species_columnreference="species",
                    source_columnreference="source",
                    site_columnreference="site",
                    normalisecomponents::Bool=true,
                    return_sites::Bool = true,
                    component_delimiter = "~|~"
                    )

    userlocations = normalize_userlocations(userlocations)
    asymmetricparams = normalize_userlocations(asymmetricparams)
    ignore_missing_singleparams = String.(ignore_missing_singleparams)
    ignore_headers = String.(ignore_headers)

    species_columnreference = String(species_columnreference)
    source_columnreference = String(source_columnreference)
    site_columnreference = String(site_columnreference)

    options = ParamOptions(;userlocations,
                            asymmetricparams,
                            ignore_missing_singleparams,
                            ignore_headers,
                            verbose,
                            species_columnreference,
                            source_columnreference,
                            site_columnreference,
                            normalisecomponents,
                            return_sites,
                            component_delimiter)

    # locations is a list of paths relative to the Clapeyron database directory.
    # userlocations is a list of paths input by the user.
    # If parameters exist ∈ multiple files, Clapeyron gives priority to files ∈ later paths.
    # asymmetricparams is a list of parameters for which matrix reflection is disabled.
    # ignore_missingsingleparams gives users the option to disable component existence check ∈ single params.
    return getparams(format_components(components),locations,options)
end

function getparams(components,locations,options::ParamOptions)
    #generate one string of params
    filepaths = flattenfilepaths(locations,options.userlocations)
    #merge all found params
    allparams,allnotfoundparams = createparams(components, filepaths, options)

    #generate sites, if any
    sites = buildsites(components,allparams,allnotfoundparams,options)# Union{SiteParam,Nothing}

    #generate ClapeyronParams
    result = compile_params(components,allparams,allnotfoundparams,sites,options)
    #check values
    for v ∈ values(result)
        is_valid_param(v,options)
    end

    return result
end

# Modified findsites: returns (sites_per_comp, any_sites, all_site_names_set)
function findsites(data::Dict, components; verbose = false)
    nc = length(components)
    sites = Dict(components .=> [Set{String}() for _ in 1:nc])
    all_sites_set = Set{String}()
    any_sites = false

    for raw in values(data)
        if raw.type === assocdata
            for (c1, c2, s1, s2) in raw.component_info
                push!(sites[c1], s1)
                push!(sites[c2], s2)
                push!(all_sites_set, s1)
                push!(all_sites_set, s2)
                any_sites = true
            end
        end
    end

    output = [sort!(collect(sites[comp])) for comp in components]
    verbose && @info("Found sites for $components are $(output).")
    return output, any_sites, all_sites_set
end


# Optimized buildsites
function buildsites(components, allparams, allnotfoundparams, options)
    options.return_sites || return nothing

    # Check if any association data exists (found or not found)
    assoc_data_found = any(x -> x.type == assocdata, values(allparams))
    assoc_data_notfound = any(x -> x == assocdata, values(allnotfoundparams))
    !assoc_data_found && !assoc_data_notfound && return nothing

    # Single pass: get sites and any_sites flag
    allcomponentsites, any_sites, all_sites_set = findsites(allparams, components)
    if !any_sites
        return SiteParam(components)   # association data exists but no actual sites
    end

    # Unique site names (sorted for deterministic order)
    v = sort(collect(all_sites_set))

    # Build mapping site => column name for number of sites
    n_sites_columns = isempty(options.n_sites_columns) ?
        Dict{String,String}(vi => string("n_", vi) for vi in v) :
        options.n_sites_columns

    nc = length(components)
    n_sites_dict = Dict{String,Vector{Int}}()
    for vi in v
        ki = n_sites_columns[vi]
        if haskey(allparams, ki)
            n_sites_dict[vi] = compile_single_vec(components, allparams[ki])
        else
            options.verbose && @warn("no columns found containing number of sites of type $(error_color(vi)). supposing zero sites")
            n_sites_dict[vi] = zeros(Int, nc)
        end
    end

    n_sites = packed_zeros!(Int,map(length,allcomponentsites))
    for i in 1:nc
        ni = n_sites[i]
        si = allcomponentsites[i]
        for a in eachindex(ni)
            sia = si[a]
            ni[a] = n_sites_dict[sia][i]
        end
    end

    # Collect source CSV paths from the single‑parameter columns used
    sourcecsvs = String[]
    for vi in v
        ki = n_sites_columns[vi]
        if haskey(allparams, ki)
            csv_vi = allparams[ki].csv
            csv_vi !== nothing && append!(sourcecsvs, csv_vi)
            delete!(allparams, ki)   # consume the single params
        end
    end
    unique!(sourcecsvs)
    return SiteParam(components, allcomponentsites, n_sites, sourcecsvs)
end

function getparams(groups::GroupParameter, locations=String[],options::ParamOptions=DefaultOptions)
    return getparams(groups.flattenedgroups, locations, options)
end



#hooks to transform arbitrary data formats into namedtuples or dicts
to_nt(x) = x

#hook to check if a struct can be transformed into a named tuple
can_nt(x) = false
can_nt(x::AbstractDict) = true
can_nt(x::NamedTuple) = true

@nospecialize
function createparams(components,
                    filepaths,
                    options::ParamOptions = DefaultOptions,
                    parsegroups = :off)

    allparams = Dict{String,RawParam}()
    allnotfoundparams = Dict{String,CSVType}()
    #in case of NamedTuple or Dict user-provided params, the filepath string should be empty.
    #but if its not, parse those anyway.

    if isempty(filepaths) && options.verbose
        @info "No string filepaths in the input."
    end


    normalised_components = normalisestring.(components,options.normalisecomponents)
    components_dict = Dict(v => k for (k,v) ∈ pairs(normalised_components))
    components_data = (components,components_dict,normalised_components)

    for raw_filepath ∈ filepaths

        _replace = startswith(raw_filepath,"@REPLACE")
        if _replace
            filepath = chop(raw_filepath,head = 9, tail = 0)
        else
            filepath = raw_filepath
        end
        csv_options = read_csv_options(filepath)
        csvtype = csv_options.csvtype

        if csvtype == groupdata && parsegroups != :group
            continue
        end

        if csvtype == assocdata && !options.return_sites
            if options.verbose
                __verbose_findparams_skipassoc(filepath)
            end
            continue
        end

        if csvtype == invaliddata
            if options.verbose
                __verbose_findparams_invaliddata(filepath)
            end
            continue
        end

        foundparams, notfoundparams = findparamsincsv(components_data,filepath,options,parsegroups,csv_options)
        merge_allparams!(allparams,allnotfoundparams,foundparams,notfoundparams,_replace)
    end

    if can_nt(options.userlocations)
        foundparams, notfoundparams = findparamsinnt(components,options,parsegroups,NT_CSV_OPTIONS)
        merge_allparams!(allparams,allnotfoundparams,foundparams,notfoundparams,false)
    end

    #clean not found params.
    for (kk,vv) ∈ allparams
        delete!(allnotfoundparams,kk)
    end

    return allparams,allnotfoundparams
end

#helper function, merges params into the main list
function merge_allparams!(allparams,allnotfoundparams,foundparams,notfoundparams,_replace)
    for vv ∈ foundparams
        kk = vv.name
        #we merge if the filepath is not set to replace the current values
        if haskey(allparams,kk) && !_replace
            vv2 = allparams[kk]
            vvx = joindata!(vv2,vv)
        else
            vvx = vv
        end
        allparams[kk] = vvx
    end
    #Merge not found data
    for (kk,vv) ∈ pairs(notfoundparams)
        if haskey(allnotfoundparams,kk)
            vv2 = allnotfoundparams[kk]
            vvx, success = joindata!(vv2,vv)
            !success && error_clashing_headers(vv2,vvx,kk) #Clashing headers error
        else
            vvx = vv
        end
        allnotfoundparams[kk] = vv
    end

    if _replace #if the parameter is not found, that means that we want to erase that param.
        for (kk,vv) ∈ pairs(notfoundparams)
            delete!(allparams,kk)
        end
    end
    return nothing
end

@specialize
function compile_params(components,allparams,allnotfoundparams,sites,options)

    #Compile Params
    result = Dict{String,ClapeyronParam}()
    for (k,v) ∈ allparams
        if !(v.type == assocdata && !options.return_sites)
            result[k] = compile_param(components,k,v,sites,options)
        end
    end
    for (kk,vv) ∈ allnotfoundparams
        result[kk] = compile_param(components,kk,vv,sites,options)
    end

    #add missing single params, if not in the input databases.
    for prop in options.ignore_missing_singleparams
        get!(result,prop) do
            options.verbose && __verbose_missing_singleparams_added(prop)
            if prop in options.asymmetricparams
                PairParam(prop,components)
            else
                SingleParam(prop,components)
            end
        end
    end

    if sites !== nothing && options.return_sites
        haskey(result,"sites") && throw(error("cannot overwrite \"sites\" key, already exists!"))
        result["sites"] = sites
    end

    return result
end

@noinline function __verbose_missing_singleparams_added(prop)
    color_prop = info_color(prop)
    @info "Optional Single Parameter $color_prop not found in any databases. Adding empty placeholder"
end

@noinline function _col_indices_error(header)
    throw(error("Header ", header, " not found."))
end

function valid_headerparams_indices(csvheaders, options::ParamOptions = DefaultOptions)
    ignorelist = deepcopy(options.ignore_headers)
    push!(ignorelist,options.species_columnreference)
    push!(ignorelist,options.source_columnreference)
    push!(ignorelist,options.site_columnreference)
    map!(normalisestring,ignorelist,ignorelist)

    v = Int[]
    for i in eachindex(csvheaders)
        header = csvheaders[i]
        filter_header = normalisestring(header; tofilter=r"[ \-\_\d]") ∉ ignorelist
        if filter_header
            push!(v,i)
        end
    end
    return v
end

function col_indices(csvtype,headernames,options=DefaultOptions)

    species = normalisestring(options.species_columnreference)
    source = normalisestring(options.source_columnreference)

    if csvtype === singledata || csvtype == groupdata
        comp1 = species
        comp2 = species
        site1 = species
        site2 = species
    elseif csvtype === pairdata || csvtype == assocdata
        comp1 = species * "1"
        comp2 = species * "2"
        if csvtype == pairdata
            site1 = species
            site2 = species
        else
            site = normalisestring(options.site_columnreference)
            site1 = site * "1"
            site2 = site * "2"
        end
    end

    idx_species = 0
    idx_species1 = 0
    idx_species2 = 0
    idx_sites1 = 0
    idx_sites2 = 0
    idx_source = 0

    for i in eachindex(headernames)
        
        header = headernames[i]

        if (csvtype === singledata || csvtype == groupdata) && idx_species == 0 && headernames[i] == species
                idx_species = i
                continue
            end
        
        if idx_species1 == 0 && header == comp1
            idx_species1 = i
            continue
        end

        if idx_species2 == 0 && header == comp2
            idx_species2 = i
            continue
        end 

        if csvtype == assocdata && idx_sites1 == 0 && header == site1
            idx_sites1 = i
            continue
        end 

        if csvtype == assocdata && idx_sites2 == 0 && header == site2
            idx_sites2 = i
            continue
        end 

        if header == source
            idx_source = i
        end
    end


    if csvtype === singledata || csvtype == groupdata
        iszero(idx_species) && _col_indices_error(species)
    elseif csvtype === pairdata || csvtype == assocdata
        iszero(idx_species1) && _col_indices_error(comp1)
        iszero(idx_species2) && _col_indices_error(comp2)
        if csvtype == assocdata
            iszero(idx_sites1) && _col_indices_error(site1)
            iszero(idx_sites2) && _col_indices_error(site2)
        end
    end

    _single = idx_species
    _pair = (idx_species1,idx_species2)
    _assoc = (idx_sites1,idx_sites2)
    return (_single,_pair,_assoc,idx_source)
end

function read_csv(filepath,options::ParamOptions,sep = :auto)::CSV.File
    #actual reading
    ignorelist = deepcopy(options.ignore_headers)
    map!(normalisestring,ignorelist,ignorelist)
    push!(ignorelist,"column") #autogenerated name by CSV.jl
    function _drop(i,name)
        norm_header = normalisestring(string(name))
        normalisestring(norm_header; tofilter=r"[ \-\_\d]") ∈ ignorelist
    end
    if sep == :auto
        sep = read_csv_options(filepath)[:sep]
    end

    _delims = (comma = ',',space = ' ')
    if sep isa Symbol
        _delim = get(_delims,sep,string(sep))
    else
        _delim = sep
    end
    if is_inline_csv(filepath)
        df = CSV.File(IOBuffer(filepath); header=3, pool=0,normalizenames=true, silencewarnings=true,drop = _drop, stringtype = String, delim = _delim, ntasks  = 1,buffer_in_memory = true)
    else
        df = CSV.File(filepath; header=3, pool=0,silencewarnings=true,drop = _drop, stringtype = String,delim = _delim, ntasks  = 1)
    end
    return df
end

function read_csv(path;relativetodatabase = true)
    path_norm = getpath(path;relativetodatabase)
    return read_csv(path_norm,ParamOptions(ignore_headers = String[]))
end

function indices_in_csv(species,comp_dict,norm,sep)
    matched_rows = Int[]
    comp_indices = Int[]
    for (row, s1) in enumerate(species)
        s1m = normalisestring(s1,norm) #normalized species
        for each_s1m in eachsplit(s1m,sep,keepempty = false) #normalized species in each entry separated by sep
            i = get(comp_dict, each_s1m, 0)
            if i != 0
                push!(matched_rows,row)
                push!(comp_indices,i)
            end
        end
    end
    return matched_rows,comp_indices,comp_indices
end

function indices_in_csv(species1,species2,comp_dict,norm,sep)
    matched_rows = Int[]
    comp1_indices = Int[]
    comp2_indices = Int[]
    for (row, s1) in enumerate(species1)
        s1m = normalisestring(s1,norm) #normalized species
        for each_s1m in eachsplit(s1m,sep,keepempty = false) #normalized species in each entry separated by sep
            i = get(comp_dict, each_s1m, 0)
            i == 0 && continue
            s2m = normalisestring(species2[row],norm)
            for each_s2m in eachsplit(s2m,sep,keepempty = false)
                j = get(comp_dict, each_s2m, 0)
                j == 0 && continue
                push!(matched_rows,row)
                push!(comp1_indices,i)
                push!(comp2_indices,j)
            end
        end
    end
    return matched_rows,comp1_indices,comp2_indices
end

function findparamsincsv(_components,filepath,
    options::ParamOptions = DefaultOptions,
    parsegroups = :off,
    csv_file_options = read_csv_options(filepath) #we do a preliminar reading of the CSV here
    )

    if _components isa Vector{String}
        normalised_components = normalisestring.(_components,options.normalisecomponents)
        components_dict = Dict(v => k for (k,v) ∈ pairs(normalised_components))
        components = _components
    else
        components,components_dict,normalised_components = _components
    end

    verbose = options.verbose
    normalisecomponents = options.normalisecomponents
    component_delimiter = options.component_delimiter
    csvtype = csv_file_options.csvtype
    no_parsegroups = parsegroups == :off
    correct_group = (parsegroups == :group && csvtype == groupdata)
    grouptype = csv_file_options.grouptype

    sep = get(csv_file_options,:sep,:comma)
    df = read_csv(filepath,options,sep)

    csvheaders = String.(Tables.columnnames(df))
    normalised_csvheaders = normalisestring.(csvheaders)
    headerparams_indices = valid_headerparams_indices(normalised_csvheaders,options) #removes all ignored header params
    headerparams = view(csvheaders,headerparams_indices)

    verbose && __verbose_findparams_start(filepath,components,headerparams,parsegroups,csvtype,grouptype)

    foundvalues = Vector{RawParam}(undef,0)
    notfoundvalues = Dict{String,CSVType}()

    csvtype == groupdata && no_parsegroups && return foundvalues, notfoundvalues
    iszero(length(headerparams)) && return foundvalues, notfoundvalues

    single_idx,pair_idx,assoc_idx,source_idx = col_indices(csvtype,normalised_csvheaders,options)
    lookupcolumnindex = single_idx
    lookupcolumnindex1,lookupcolumnindex2 = pair_idx
    lookupsitecolumnindex1,lookupsitecolumnindex2 = assoc_idx

    lookupcolumnindex = max(lookupcolumnindex,lookupcolumnindex1)

    if csvtype == singledata || correct_group
        species1_list   = Tables.getcolumn(df,lookupcolumnindex)
        species2_list   = species1_list
        site1_list      = species1_list
        site2_list      = species1_list
        matched_rows,comp1_indices,comp2_indices = indices_in_csv(species1_list,components_dict,normalisecomponents,component_delimiter)
    elseif csvtype == pairdata || csvtype == assocdata
        species1_list   = Tables.getcolumn(df,lookupcolumnindex1)
        species2_list   = Tables.getcolumn(df,lookupcolumnindex2)
        site1_list      = csvtype == pairdata ? species1_list : Tables.getcolumn(df,lookupsitecolumnindex1)
        site2_list      = csvtype == pairdata ? species1_list : Tables.getcolumn(df,lookupsitecolumnindex2)
        matched_rows,comp1_indices,comp2_indices = indices_in_csv(species1_list,species2_list,components_dict,normalisecomponents,component_delimiter)
    else
        species1_list   = Tables.getcolumn(df,lookupcolumnindex)
        species2_list   = species1_list
        site1_list      = species1_list
        site2_list      = species1_list
        matched_rows    = Int[]
        comp1_indices,comp2_indices = matched_rows,matched_rows
    end

    N_found = length(matched_rows)

    EMPTY_STR = ""

    if source_idx != 0 && N_found > 0
        all_sources = Vector{String}(undef,N_found)
        source_col = Tables.getcolumn(df,source_idx)
        map!(Base.Fix2(coalesce,EMPTY_STR),all_sources,view(source_col,matched_rows))
    else
        all_sources = String[]
    end

    for (i,icol) in enumerate(headerparams_indices)
        header = headerparams[i]
        if N_found == 0
            notfoundvalues[strip(header)] = csvtype
            continue
        end
    
        col = Tables.getcolumn(df, icol)
        values_raw = col[matched_rows]
        N_nonmissing = count(!ismissing,values_raw)

        if N_nonmissing == 0
            notfoundvalues[strip(header)] = csvtype
            continue
        end

        values = nonmissingvec(values_raw)

        header          = headerparams[i]
        component_info  = Vector{NTuple{4,String}}(undef,N_nonmissing)
        sources         = Vector{String}(undef,N_nonmissing)
        sourcecsvs      = Vector{String}(undef,N_nonmissing)
        sourcecsvs      .= filepath
        j = 0
        for _j in 1:N_found
            vj = values_raw[_j]
            ismissing(vj) && continue
            j += 1
            jx = matched_rows[_j]
            if csvtype == singledata || correct_group
                c = comp1_indices[_j]
                component_info[j] = (components[c],components[c],EMPTY_STR,EMPTY_STR)
            elseif csvtype == pairdata && no_parsegroups
                c1,c2 = comp1_indices[_j],comp2_indices[_j]
                component_info[j] = (components[c1],components[c2],EMPTY_STR,EMPTY_STR)
            elseif csvtype == assocdata && no_parsegroups
                c1,c2 = comp1_indices[_j],comp2_indices[_j]
                component_info[j] = (components[c1],components[c2],site1_list[jx],site2_list[jx])
            end
            sources[j] = length(all_sources) != 0 ? all_sources[_j] : EMPTY_STR
        end
        raw = RawParam(string(strip(header)),component_info,values,sources,sourcecsvs,csvtype,grouptype)

        push!(foundvalues, raw)
    end

    verbose && __verbose_findparams_found(foundvalues) #print all found values
    verbose && __verbose_findparams_not_found(notfoundvalues) #print all found values

    return foundvalues, notfoundvalues
end


#find params in named tuple, transforms form named tuple to Dict{RawParam}
function findparamsinnt(components,
    options::ParamOptions,
    parsegroups = :off,
    csv_file_options = NT_CSV_OPTIONS) #default options

    nt = to_nt(options.userlocations)
    foundvalues = Vector{RawParam}(undef,0)
    notfoundvalues = Dict{String,CSVType}()
    #this algorithm is less strict that what we have in CSVs. but allows us to parse named tuples

    for (k,v) in pairs(nt)
        ks = string(k)
        if ks == "groups" && parsegroups == :groups
            param = RawParam(ks,nothing,copy(v),nothing,nothing,groupdata,:unknown)
            push!(foundvalues,param)
        elseif (ks == "epsilon_assoc" || ks == "bondvol") && parsegroups == :off && v === nothing #TODO: what to do here in case of other assoc names?
            notfoundvalues[ks] = assocdata
        elseif v isa AbstractVector && parsegroups == :off
            vv = convert(Vector,v)
            param = RawParam(ks,nothing,copy(vv),nothing,nothing,singledata,:unknown)
            push!(foundvalues,param)
        elseif v isa AbstractMatrix && parsegroups == :off
            vv = vec(convert(Matrix,v))
            param = RawParam(ks,nothing,copy(vv),nothing,nothing,pairdata,:unknown)
            push!(foundvalues,param)
        elseif v isa Number && parsegroups == :off && length(components) == 1
            param = RawParam(ks,nothing,[v],nothing,nothing,singledata,:unknown)
            push!(foundvalues,param)
        elseif v isa AbstractDict
            val1 = first(values(v))
            assoc_values = Vector{typeof(val1)}(undef,0)
            param = RawParam(ks,Vector{NTuple{4,String}}(undef,0),assoc_values,String[],String[],assocdata,:unknown)
            empty_string = ""
            for (k_dict,v_dict) in pairs(v)
                sp1,s1 = first(k_dict)
                sp2,s2 = last(k_dict)
                push!(param.component_info,(sp1,sp2,s1,s2))
                push!(param.data,v_dict)
                push!(param.sources,empty_string)
                push!(param.csv,empty_string)
                push!(foundvalues,param)
            end
        else
            throw(error("cannot parse combination key = $k, value = $v as a valid parameter."))
        end
    end

   # verbose && __verbose_findparams_found(foundvalues) #print all found values

    return foundvalues, notfoundvalues
end

function nonmissingvec(X::AbstractVector{T}) where T
    Y = Vector{nonmissingtype(T)}(undef,length(X))
    k = 0
    for i in eachindex(X)
        if !ismissing(X[i])
            k += 1
            Y[k] = X[i]
        end
    end
    resize!(Y,k)
    return Y
end

#verbose functionality, is executed for each csv when verbose == true

function __verbose_findparams_invaliddata(filepath)
    @warn "Skipping $filepath, cannot infer correct csv type. Check line 2 of the CSV to see if it has valid information."
end

function __verbose_findparams_skipassoc(filepath)
    @warn "Skipping association file $filepath, the option return_sites was set to false."
end


function __assoc_string(pair)
    "($(pair[1]),$(pair[3])) >=< ($(pair[2]), $(pair[4]))"
end

function __verbose_findparams_start(filepath,components,headerparams,parsegroups,csvtype,grouptype)
    csv_string = Symbol(csvtype)
    no_parsegroups = parsegroups == :off
    if no_parsegroups
        if csvtype == groupdata
            @info("Skipping $csv_string csv $filepath")
        elseif iszero(length(headerparams))
            @info("Skipping $csv_string csv at $filepath, couldn't find any valid headers ...")
        else
            @info("Searching for $csv_string headers $headerparams for query $components at $filepath ...")
        end
    else
        if csvtype == groupdata
            @info("Searching for groups for components $components at $filepath ...")
        else
            @info("Skipping $csv_string csv at $filepath")
        end
    end
    if grouptype != :unknown
        @info("group type: $grouptype")
    end
end

function __verbose_findparams_not_found(notfoundvalues)
    for k ∈ keys(notfoundvalues)
        kk = info_color(k)
        @info("No data found for $kk")
    end
end

function __verbose_findparams_found(foundvalues)
    for v ∈ foundvalues
        if v.type == singledata
            io = IOBuffer()
            show_pairs(io,first.(v.component_info),v.data," => ",quote_string = false)
            vals = String(take!(io))
            kk = info_color(v.name)
            TT = eltype(v.data)
            @info("""Found single component data: $kk with $TT values:
            $vals
            """)
        elseif v.type == pairdata
            io = IOBuffer()
            show_pairs(io,first.(v.component_info,2),v.data," => ",quote_string = false)
            vals = String(take!(io))
            kk = info_color(v.name)
            TT = eltype(v.data)
            @info("""Found pair component data: $kk with $TT values:
            $vals
            """)
        elseif v.type == assocdata
            io = IOBuffer()
            show_pairs(io,__assoc_string.(v.component_info),v.data," => ",quote_string = false)
            vals = String(take!(io))
            kk = info_color(v.name)
            TT = eltype(v.data)
            @info("""Found association component data: $kk with $TT values:
            $vals
            """)
        elseif v.type == groupdata
            #println(val) for val ∈ v.data)
            vals = "Dict("
            for (pair,val) ∈ zip(v.component_info,v.data)
                pairi = pair[1] * " => " * string(val)
                vals = vals * pairi * ", "
            end
            vals = chop(vals,tail=2) *")"
            @info("""Found group data:
            $vals
            """)
            #@info("TODO: parse intragroup data for debug")
        end
    end
end

const READCSVTYPE_KEYWORDS  = Set(["like", "single", "unlike", "pair", "assoc", "association", "group", "groups","intragroup","intragroups"])


function read_csv_options(filepath::AbstractString)
    return _read_csv_options(getline(String(filepath), 2))
end
#=
function read_csv_options(filepath)
    return 2
end=#

function _read_csv_options(line::String)
    vec_re = r"\[.*\]"
    maybe_opts_vec = match(vec_re,line)
    json_re = r"\{.*\}"
    maybe_opts_json = match(json_re,line)

    # Searches for type from second line of CSV.
    has_csv_options_vec = !isnothing(maybe_opts_vec)
    has_csv_options_json = !isnothing(maybe_opts_json)
    if has_csv_options_json
        __get_options(maybe_opts_json.match,:json)
    elseif has_csv_options_vec
        opts = chop(maybe_opts_vec.match,head = 1,tail = 1)
        return __get_options(opts,:vec)
    else
        data = [""]
        words = split(lowercase(strip(line, ',')), ' ')

        maybe_csvdata = false
        for word in words
            if word in READCSVTYPE_KEYWORDS && !maybe_csvdata
                maybe_csvdata = true
                data[1] = word
            elseif word in READCSVTYPE_KEYWORDS && maybe_csvdata
                data[1] = ""
            end
        end
        return (csvtype = _readcsvtype(data[1]),grouptype = :unknown, sep = :comma)
    end
end

const NT_CSV_OPTIONS = (csvtype = namedtupledata,grouptype = :unknown,sep = :comma)

function _readcsvtype(collection)
    length(collection) != 1 && return invaliddata
    key = only(collection)
    return _readcsvtype(key)
end

function _readcsvtype(key::AbstractString)
    key == "single" && return singledata
    key == "like" && return singledata
    key == "pair" && return pairdata
    key == "unlike" && return pairdata
    key == "assoc" && return assocdata
    key == "group" && return groupdata
    key == "groups" && return groupdata
    key == "intragroup" && return groupdata
    key == "intragroups" && return groupdata
    key == "invalid" && return invaliddata
    return invaliddata
end

"""
    parse_bracket_format(s) -> Dict{SubString{String}, Tuple{Bool,SubString{String}}}

Parse a bracket-encoded key-value string. Used in Clapeyron.jl CSV file headers.
returns a dictionary of key-value pairs. each value is a tuple of a boolean that indicates if the text is a vector of values
and the text content. parsing of vectors is optional.
"""
function parse_bracket_format(s::AbstractString,brackets_removed = false)
    if !brackets_removed
        s = strip(s)
        m = match(r"^\[(.*)\]$"s, s)          # regex: strip outer []
        m === nothing && throw(ArgumentError("Input must be enclosed in '[…]': $(repr(s))"))
        content = strip(m[1])
    else
        content = s
    end
    isempty(content) && return Dict{String, Union{String, Vector{String}}}()


    result = Dict{SubString{String}, Tuple{Bool,SubString{String}}}()
    seg_start = firstindex(s)
    done = false
    while !done
        add,idx,seg_start,done = _split_pairs(content,seg_start)
        if add
            ss = strip(SubString(content,idx))
            k,v = _parse_pair(ss)
            result[k] = v
        end
    end
    return result
end

"""
    _split_pairs(s) -> Vector{String}

Split the flat content string into individual "key=value" segments, considering quotes and whitespace.

"""
function _split_pairs(s::AbstractString,seg_start)
    in_q  = false
    depth = 0
    for i in seg_start:lastindex(s)
        c = s[i]
        if in_q
            c == '"' && (in_q = false)
        else
            if     c == '"';                in_q  = true
            elseif c == '(';               depth += 1
            elseif c == ')';               depth -= 1
            elseif c == ',' && depth == 0
                indices0 = seg_start:prevind(s, i)
                seg = strip(SubString(s,indices0))
                if !all(isspace,seg)
                    return true,indices0,nextind(s, i),false
                end
            end
        end
    end

    indices_end = seg_start:lastindex(s)
    tail = strip(SubString(s,indices_end))
    if !all(isspace,tail)
        return true,indices_end,0,true
    else
        return false,indices_end,0,true
    end
end

"""
    _parse_pair(s) -> (key::SubString{String}, Tuple{Bool,SubString{String}})

Split one "key = value" segment at the first '=' that is not inside a quoted
string or a parenthesised group.
"""
function _parse_pair(s::AbstractString)
    in_q  = false
    depth = 0
    eq_i  = nothing

    for i in eachindex(s)
        c = s[i]
        if in_q
            c == '"' && (in_q = false)
        else
            if     c == '"';                in_q  = true
            elseif c == '(';               depth += 1
            elseif c == ')';               depth -= 1
            elseif c == '=' && depth == 0; eq_i = i; break
            end
        end
    end

    eq_i === nothing && throw(ArgumentError("No '=' found in segment: $(repr(s))"))

    key = strip(SubString(s,firstindex(s):prevind(s, eq_i)))
    val = strip(SubString(s,nextind(s, eq_i):lastindex(s)))
    count_quotes = count(isequal('\"'),val)
    if startswith(val,'\"') && endswith(val,'\"') && count_quotes == 2
        val = strip(isequal('\"'),val)
    end
    
    is_vec = count_quotes > 2 || (startswith(val, '(') && endswith(val, ')'))
    return key, (is_vec,val)
end

"""
    _parse_vec(s) -> Vector{SubString{String}}

Convert a string into a vector of strings acording to the following rules:
- Each quoted value is an individual value, if there are more than 2 quotes present.
- Text between parentheses is always a vector with one or more quoted strings
- If no quotes and no parentheses are found, each word separated by whitespace is an individual value.
"""
function _parse_vec(s::AbstractString)
    reg = r"\"([^\"]*)\""
    # ── parenthesised list ─────────────────────────────────────────────────
    if startswith(s, '(') && endswith(s, ')')
        inner = SubString(s,nextind(s, firstindex(s)):prevind(s, lastindex(s)))
        return [m[1] for m in eachmatch(reg, inner)]
    end
    # ── one or more quoted tokens ──────────────────────────────────────────
    if count(isequal('\"'),s) > 2
        return [m[1] for m in eachmatch(reg, s)]
    end

    #split by spaces
    return split(s)
end

function __get_options(data,type)
    if type == :vec
        opts_dict = parse_bracket_format(data,true)
        a,b,c = "csvtype","grouptype","sep"
        _csvtype = if haskey(opts_dict,a)
            _readcsvtype(opts_dict[SubString(a)][2])
        else
            invaliddata
        end

        _grouptype = if haskey(opts_dict,b)
            Symbol(opts_dict[SubString(b)][2])
        else
            :unknown
        end

        _sep = if haskey(opts_dict,c)
            Symbol(opts_dict[SubString(c)][2])
        else
            :comma
        end
        return (csvtype = _csvtype,grouptype = _grouptype,sep = _sep)
    elseif type == :json
        json_dict = JSON.parse(data)
        _csvtype = _readcsvtype(get(json_dict,"csvtype","invalid"))
        _grouptype = Symbol(get(json_dict,"grouptype","unknown"))
        #_estimator = Symbol(get(json_dict,"method","error"))
        #maybe_species = get(json_dict,"species","all")
        #_species = maybe_species isa AbstractString ? [String(maybe_species)] : String.(maybe_species)
        _sep = Symbol(get(json_dict,"sep","comma"))
        return (csvtype = _csvtype,grouptype = _grouptype, sep = _sep)
    else
        throw(error("Clapeyron.__get_options: invalid type. expected :json or :vec, got $type"))
    end
end

include("database_group.jl")
include("database_writer.jl")

export getparams
