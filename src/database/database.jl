@enum CSVType invaliddata namedtupledata singledata pairdata assocdata groupdata structgroupdata
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
    params, sites = getparams(components,locations;kwargs...)
returns a `Dict{String,ClapeyronParam}` containing all the parameters found for the list of components
in the available CSVs. `locations` are the locations relative to `Clapeyron` database. the available keywords are the ones used ∈ [`ParamOptions`](@ref)
if `return_sites` is set to false, `getparams` will only return the found params.

## Single to Pair promotion

When reading multiple CSVs, if a parameter name appears in a single paramter file and in a pair parameter file, the single parameter values will be promoted to be the diagonal values of the pair interaction matrix:

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
This promotion fails only happens in Single-Pair combinations. it fails otherwise.

## In-memory CSV parsing

if you pass any string starting with `Clapeyron Database File`, it will be parsed as a CSV instead of being used as a filepath:

```julia
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
In a way, is a path shortcut used internally by Clapeyron to parse it's own database. you can change the path where `@DB` points to (or add other path shortcuts), via adding a corresponding entry to the `Clapeyron.SHORT_PATHS` Dict.
- `@REPLACE`: Any filepath starting with `@REPLACE` will clear all previous appearances of the parameter names found in the CSV that contains the prefix.
- `@REMOVEDEFAULTS`: it is used alone, and needs to be passed at the first position of the vector of `userlocations`. it will skip parsing of the default parameters:

The effect of the the parser can be summarized by the following examples:

```
model = PCSAFT(["water"],userlocations = ["@REMOVEDEFAULTS"]) #fails, no parameters found, no CSV parsed
model = PCSAFT(["water"],userlocations = ["@REPLACE/empty_params.csv"]) #fails, no parameters found, default parameters parsed and then removed
model = PCSAFT(["water"],userlocations = ["@REPLACE/my_pcsaft_kij.csv"]) #success, default kij parameters replaced by the ones on `my_pcsaft_kij.csv`
model = PCSAFT(["water"],userlocations = ["@REMOVEDEFAULTS","@DB/SAFT/PCSAFT","@DB/properties/molarmass.csv"]) #sucess. default parameters csv removed, and parsed again, using the @DB prefix to point to the default database.
```

You can use the `@REPLACE` keyword in a in-memory CSV by adding it at the start of the string, followed by an space:
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

The second line of the csv is used for comments and to identify the type of CSV used. for example:

    ```
x = \"\"\"Clapeyron Database File
       in memory like parameters
       species,a,b
       sp1,1000,0.05
       sp2,700,0.41
       \"\"\"
```
Will be parsed as a table with single parameter data. if you want more flexibility, you can instead pass the csvtype between brackets:

```
x = \"\"\"Clapeyron Database File
       i can write anything here, unlike, association [csvtype = like] but the csv type is already specified.
       species,a,b
       sp1,1000,0.05
       sp2,700,0.41
       \"\"\"
```
additionaly, there are some cases when you want to absolutely sure that your types don't clash with the default values. this is the case with different group parametrizations of UNIFAC (Dormund, VTPR, PSRK):

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

The models are the same (`UNIFAC`), but the group parametrizations are different. this is specified with the `grouptype` keyword. for example, if we see `UNIFAC_groups.csv`, it starts with:

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
       paramterization 1 [csvtype = like,grouptype = param1]
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

If we pass the same parameters, with different group types, the parser will fail

```julia-repl
julia> Clapeyron.getparams(["sp1","sp2"],userlocations = [x1,x2])
ERROR: cannot join two databases with different group types:
current group type: param1
incoming group type: fitted
```

Note, that the parser will not fail if you pass different parameters with different group types (For example if `a` has `param1` group type and `b` has `fit` group type)
"""
function getparams(components,
                    locations::Array{String,1}=String[];
                    userlocations = String[],
                    asymmetricparams::Vector{String}=String[],
                    ignore_missing_singleparams::Vector{String}=String[],
                    ignore_headers::Vector{String} = IGNORE_HEADERS,
                    verbose::Bool=false,
                    species_columnreference::String="species",
                    source_columnreference::String="source",
                    site_columnreference::String="site",
                    group_columnreference::String="groups",
                    normalisecomponents::Bool=true,
                    return_sites::Bool = true,
                    component_delimiter::String = "~|~"
                    )

    options = ParamOptions(;userlocations,
                            asymmetricparams,
                            ignore_missing_singleparams,
                            ignore_headers,
                            verbose,
                            species_columnreference,
                            source_columnreference,
                            site_columnreference,
                            group_columnreference,
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

function getparams(components::Vector{String},locations::Vector{String},options::ParamOptions)
    #generate one string of params
    filepaths = flattenfilepaths(locations,options.userlocations)
    #merge all found params
    allparams,allnotfoundparams = createparams(components, filepaths, options)
    #generate ClapeyronParams
    result, allcomponentsites = compile_params(components,allparams,allnotfoundparams,options)
    #check values
    for v ∈ values(result)
        is_valid_param(v,options)
    end

    if !(options.return_sites)
        return result
    end
    if any(x isa AssocParam for x ∈ values(result))
        sites = buildsites(result,components,allcomponentsites,options)
        if !haskey(result,"sites")
            result["sites"] = sites
            return result
        else
            throw(error("cannot overwrite \"sites\" key, already exists!"))
        end
    else
        return result
    end
end

function buildsites(result,components,allcomponentsites,options)


    v = String[]
    for sitei ∈ allcomponentsites
        append!(v,sitei)
    end
    unique!(v)
    #no sites encountered in assoc params
    iszero(length(v)) && return SiteParam(components)

    #=
    build our own dict, by using the transformation X => n_X
    =#
    if isempty(options.n_sites_columns)
        n_sites_columns = Dict{String,String}(vi => string("n_",vi) for vi in v)
    else
        n_sites_columns = options.n_sites_columns
    end

    #check for missing number of sites
    if !any(haskey(result,get(n_sites_columns,vi,"")) for vi ∈ v)
       return __warning_no_site_vals(result,components)
    end

    n_sites_dict = Dict{String,SingleParam{Int}}(vi => result[n_sites_columns[vi]] for vi ∈ v)
    return SiteParam(n_sites_dict,allcomponentsites)
end

function __warning_no_site_vals(result,components)
    @error """No columns containing number of sites were found. Supposing zero sites.
    If your model doesn't use sites, but some input CSV folders contain Association params. consider using return_sites=false.
    If there are columns containing number of sites, then those aren't recognized.
    Consider passing a Dict with the mappings between the sites and the name of the column containing the number of said sites,
    with the n_sites_columns keyword"
    """
    assoc_csv = Set(String[])
    for x ∈ values(result)
        if x isa AssocParam
            for csvx ∈ x.sourcecsvs
            push!(assoc_csv,csvx)
            end
        end
    end
    assoc_csv = collect(assoc_csv)
    @error "Parsed Association CSV were:"
    for csv ∈ assoc_csv
        println(csv)
    end
    return SiteParam(components)
end

function getparams(groups::GroupParameter, locations::Vector{String}=String[],options::ParamOptions=DefaultOptions)
    return getparams(groups.flattenedgroups, locations,options)
end

function findsites(data::Dict,components::Vector;verbose = false)
    sites = Dict(components .=> [Set{String}() for _ ∈ 1:length(components)])
    for raw ∈ values(data)
        if raw.type === assocdata
            for (c1,c2,s1,s2) ∈ raw.component_info
            push!(sites[c1], s1)
            push!(sites[c2], s2)
            end
        end
    end
    output = Array{Array{String,1}}(undef, 0)
    for component ∈ components
        push!(output, collect(sites[component]))
    end
    verbose && @info("Found sites for $components are $(output).")
    return output
end

#hooks to transform arbitrary data formats into namedtuples or dicts
to_nt(x) = x

#hook to check if a struct can be transformed into a named tuple
can_nt(x) = false
can_nt(x::AbstractDict) = true
can_nt(x::NamedTuple) = true

@nospecialize
function createparams(components::Vector{String},
                    filepaths::Vector{String},
                    options::ParamOptions = DefaultOptions,
                    parsegroups = :off)

    allparams = Dict{String,RawParam}()
    allnotfoundparams = Dict{String,CSVType}()
    #in case of NamedTuple or Dict user-provided params, the filepath string should be empty.
    #but if its not, parse those anyway.
    for filepath ∈ filepaths

        _replace = startswith(filepath,"@REPLACE")
        if _replace
            filepath = chop(filepath,head = 9, tail = 0)
        end
        csv_options = read_csv_options(filepath)
        csvtype = csv_options.csvtype

        if csvtype == groupdata && parsegroups != :group
            continue
        end

        if csvtype == structgroupdata && parsegroups != :intragroup
            continue
        end

        if csvtype == invaliddata
            if options.verbose
                __verbose_findparams_invaliddata(filepath)
            end
            continue
        end

        foundparams, notfoundparams = findparamsincsv(components,filepath,options,parsegroups,csv_options)
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
            vv = joindata!(vv2,vv)
        end
        allparams[kk] = vv
    end
    #Merge not found data
    for (kk,vv) ∈ pairs(notfoundparams)
        if haskey(allnotfoundparams,kk)
            vv2 = allnotfoundparams[kk]
            vv, success = joindata!(vv2,vv)
            !success && error_clashing_headers(vv2,vv,kk) #Clashing headers error
        end
        allnotfoundparams[kk] = vv
    end

    if _replace #if the paramter is not found, that means that we want to erase that param.
        for (kk,vv) ∈ pairs(notfoundparams)
            delete!(allparams,kk)
        end
    end
    return nothing
end
@specialize
function compile_params(components,allparams,allnotfoundparams,options)
    #Generate component sites with the RawParam tapes
    allcomponentsites = findsites(allparams,components)
    #Compile Params
    result = Dict{String,ClapeyronParam}(k => compile_param(components,k,v,allcomponentsites,options) for (k,v) ∈ allparams)
    for (kk,vv) ∈ allnotfoundparams
        result[kk] = compile_param(components,kk,vv,allcomponentsites,options)
    end
    return result,allcomponentsites
end

@noinline function _col_indices_error(header)
    throw(error("Header ", header, " not found."))
end

function col_indices(csvtype,headernames,options=DefaultOptions)
    columnreference = options.species_columnreference
    normalised_columnreference = normalisestring(columnreference)

    idx_species = 0
    idx_groups = 0
    idx_species1 = 0
    idx_species2 = 0
    idx_sites1 = 0
    idx_sites2 = 0

    if csvtype === singledata || csvtype ∈ (groupdata,structgroupdata)
        lookupcolumnindex = findfirst(isequal(normalised_columnreference), headernames)
        isnothing(lookupcolumnindex) && _col_indices_error(normalised_columnreference)
        idx_species = lookupcolumnindex
        if csvtype == groupdata
            groupcolumnreference = options.group_columnreference
            normalised_groupcolumnreference = normalisestring(groupcolumnreference)
            lookupgroupcolumnindex = findfirst(isequal(normalised_groupcolumnreference), headernames)
            isnothing(lookupgroupcolumnindex) && _col_indices_error(normalised_groupcolumnreference)
            idx_groups = lookupgroupcolumnindex
        else
            idx_groups = 0
        end

    elseif csvtype === pairdata || csvtype == assocdata
        normalised_columnreference1 = normalised_columnreference * '1'
        normalised_columnreference2 = normalised_columnreference * '2'
        lookupcolumnindex1 = findfirst(isequal(normalised_columnreference1), headernames)
        lookupcolumnindex2 = findfirst(isequal(normalised_columnreference2), headernames)
        isnothing(lookupcolumnindex1) && _col_indices_error(normalised_columnreference1)
        isnothing(lookupcolumnindex2) && _col_indices_error(normalised_columnreference2)
        idx_species1 = lookupcolumnindex1
        idx_species2 = lookupcolumnindex2
        if csvtype == assocdata
            sitecolumnreference = options.site_columnreference
            normalised_sitecolumnreference = normalisestring(sitecolumnreference)
            normalised_sitecolumnreference1 = normalised_sitecolumnreference * '1'
            normalised_sitecolumnreference2 = normalised_sitecolumnreference * '2'
            lookupsitecolumnindex1 = findfirst(isequal(normalised_sitecolumnreference1), headernames)
            lookupsitecolumnindex2 = findfirst(isequal(normalised_sitecolumnreference2), headernames)
            isnothing(lookupsitecolumnindex1) && _col_indices_error(normalised_sitecolumnreference1)
            isnothing(lookupsitecolumnindex2) && _col_indices_error(normalised_sitecolumnreference2)
            idx_sites1 = lookupsitecolumnindex1
            idx_sites2 = lookupsitecolumnindex2
        end
    end

    _single = (idx_species,idx_groups)
    _pair = (idx_species1,idx_species2)
    _assoc = (idx_sites1,idx_sites2)
    return (_single,_pair,_assoc)
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
        df = CSV.File(IOBuffer(filepath); header=3, pool=0,silencewarnings=true,drop = _drop, stringtype = String, delim = _delim, ntasks  = 1,buffer_in_memory = true)
    else
        df = CSV.File(filepath; header=3, pool=0,silencewarnings=true,drop = _drop, stringtype = String,delim = _delim, ntasks  = 1)
    end
    return df
end

function read_csv(path;relativetodatabase = true)
    path_norm = getpath(path;relativetodatabase)
    return read_csv(path_norm,ParamOptions(ignore_headers = String[]))
end

function findparamsincsv(components,filepath,
    options::ParamOptions = DefaultOptions,
    parsegroups = :off,
    csv_file_options = read_csv_options(filepath) #we do a preliminar reading of the CSV here
    )

    sourcecolumnreference = options.source_columnreference
    verbose = options.verbose
    normalisecomponents = options.normalisecomponents
    component_delimiter = options.component_delimiter
    csvtype = csv_file_options.csvtype
    no_parsegroups = parsegroups == :off
    correct_group = (parsegroups == :group && csvtype == groupdata) || (parsegroups == :intragroup && csvtype == structgroupdata)
    grouptype = csv_file_options.grouptype

    sep = get(csv_file_options,:sep,:comma)
    df = read_csv(filepath,options,sep)

    csvheaders = String.(Tables.columnnames(df))
    headerparams = valid_headerparams(csvheaders,options) #removes all ignored header params

    normalised_components = normalisestring.(components,normalisecomponents)
    components_dict = Dict(v => k for (k,v) ∈ pairs(normalised_components))

    normalised_csvheaders = normalisestring.(csvheaders)
    normalised_headerparams = normalisestring.(headerparams)

    if normalised_headerparams ⊈ normalised_csvheaders
        error("Headers ", setdiff(normalised_headerparams, normalised_csvheaders), " not present ∈ csv header.")
    end

    #function output
    foundvalues = Vector{RawParam}(undef,0)
    notfoundvalues = Dict{String,CSVType}(strip(headerparam) => csvtype for headerparam ∈ headerparams)

    normalised_sourcecolumnreference = normalisestring(sourcecolumnreference)
    getsources = false
    if normalised_sourcecolumnreference ∈ normalised_csvheaders
        getsources = true
        sourcecolumn = findfirst(isequal(normalised_sourcecolumnreference), normalised_csvheaders)::Int
    else
        sourcecolumn = 0
    end

    single_idx,pair_idx,assoc_idx = col_indices(csvtype,normalised_csvheaders,options)
    lookupcolumnindex,groupindex = single_idx
    lookupcolumnindex1,lookupcolumnindex2 = pair_idx
    lookupsitecolumnindex1,lookupsitecolumnindex2 = assoc_idx
    headerparams_indices = zeros(Int,length(normalised_headerparams))
    map!(i -> findfirst(isequal(i),normalised_csvheaders)::Int,headerparams_indices,normalised_headerparams)
    #headerparams_indices = [findfirst(isequal(i),normalised_csvheaders) for i ∈ normalised_headerparams]
    lookupcolumnindex = max(lookupcolumnindex,lookupcolumnindex1)

    verbose && __verbose_findparams_start(filepath,components,headerparams,parsegroups,csvtype,grouptype)
    #list of all species
    species_list::Vector{String} = normalisestring.(Tables.getcolumn(df,lookupcolumnindex),normalisecomponents)

    #indices where data could be (they could be missing)
    #on pair and assoc, this is just the first component, we need to reduce the valid indices again
    found_indices0,comp_indices = _indexin(components_dict,species_list,component_delimiter,1:length(species_list))
    dfR = df
    EMPTY_STR = ""
    if csvtype == singledata || correct_group
        found_indices = found_indices0
        l = length(found_indices)
        _data = dfR[found_indices]
        _sources = fill(EMPTY_STR,l)
        _csv = fill(filepath,l)
        _comp = Vector{NTuple{4,String}}(undef,l)
        for li ∈ 1:l
            _c = components[comp_indices[li]]
            _comp[li] = (_c,EMPTY_STR,EMPTY_STR,EMPTY_STR)
        end

    elseif csvtype == pairdata && no_parsegroups
        species2_list::Vector{String} = normalisestring.(Tables.getcolumn(df,lookupcolumnindex2)[found_indices0],normalisecomponents)
        found_indices2,comp_indices2 = _indexin(components_dict,species2_list,component_delimiter,1:length(species2_list))
        comp_indices1 = comp_indices[found_indices2]
        found_indices2 = found_indices0[found_indices2]
        l = length(found_indices2)
        _data = dfR[found_indices2]
        _comp = [(components[c1],components[c2],EMPTY_STR,EMPTY_STR) for (c1,c2) ∈ zip(comp_indices1,comp_indices2)]
        _sources = fill(EMPTY_STR,l)
        _csv = fill(filepath,l)
    elseif csvtype == assocdata && no_parsegroups
        species2_list = normalisestring.(Tables.getcolumn(df,lookupcolumnindex2)[found_indices0],normalisecomponents)
        found_indices2,comp_indices2 = _indexin(components_dict,species2_list,component_delimiter,1:length(species2_list))
        comp_indices1 = comp_indices[found_indices2]
        found_indices2 = found_indices0[found_indices2]
        l = length(found_indices2)
        _data = dfR[found_indices2]
        _site1 = Vector{String}(undef,l)
        _site2 = similar(_site1)
        _comp = Vector{NTuple{4,String}}(undef,l)
        for li ∈ 1:l
            _s1 = _data[li][lookupsitecolumnindex1]
            _s2 = _data[li][lookupsitecolumnindex2]
            _c1 = components[comp_indices1[li]]
            _c2 = components[comp_indices2[li]]
            _site1[li] = _s1
            _site2[li] = _s2
            _comp[li] = (_c1,_c2,_s1,_s2)
        end
        _sources = fill(EMPTY_STR,l)
        _csv = fill(filepath,l)

    elseif csvtype ∈ (groupdata,structgroupdata) && no_parsegroups
        return foundvalues, notfoundvalues
    else
        error("Filepath $filepath is of type ", string(csvtype), " and cannot be read with this function.")
    end

        #if getsources, then we actually put the sources ∈ inside the preallocated _sources vector
    if getsources
        _fill_sources!(_sources,getindex.(_data,sourcecolumn),EMPTY_STR)
    end
    #with the raw data preallocated, we now store it ∈ a RawParam.
    for (headerparam,idx) ∈ zip(headerparams,headerparams_indices)
        _vals = getindex.(_data,idx)
        raw::Any = build_raw_param(headerparam,_comp,_vals,_sources,_csv,csvtype,grouptype)
        if !iszero(length(raw))
            push!(foundvalues,raw)
        end
    end

    #store all headers that didn't had a result.
    for rawparam ∈ foundvalues
        delete!(notfoundvalues,rawparam.name)
    end

    verbose && __verbose_findparams_found(foundvalues) #print all found values

    return foundvalues, notfoundvalues
end


#find params in named tuple, transforms form named tuple to Dict{RawParam}
function findparamsinnt(components,
    options::ParamOptions,
    parsegroups = :off,
    csv_file_options = NT_CSV_OPTIONS) #default options

    verbose = options.verbose
    nt = to_nt(options.userlocations)
    foundvalues = Vector{RawParam}(undef,0)
    notfoundvalues = Dict{String,CSVType}()
    #this algorithm is less strict that what we have in CSVs. but allows us to parse named tuples

    for (k,v) in pairs(nt)
        ks = string(k)
        if k == :groups && parsegroups == :groups
            param = RawParam(ks,nothing,copy(v),nothing,nothing,groupdata,:unknown)
            push!(foundvalues,param)
        elseif k == :intragroups && parsegroups == :structgroups
            param = RawParam(ks,nothing,copy(v),nothing,nothing,structgroupdata,:unknown)
        elseif (k == :epsilon_assoc || k == :bondvol) && parsegroups == :off && v === nothing
            notfoundvalues[ks] = assocdata
        elseif v isa Vector && parsegroups == :off
            param = RawParam(ks,nothing,copy(v),nothing,nothing,singledata,:unknown)
            push!(foundvalues,param)
        elseif v isa Matrix && parsegroups == :off
            param = RawParam(ks,nothing,vec(copy(v)),nothing,nothing,pairdata,:unknown)
            push!(foundvalues,param)
        elseif v isa Number && parsegroups == :off && length(components) == 1
            param = RawParam(ks,nothing,[v],nothing,nothing,singledata,:unknown)
            push!(foundvalues,param)
        elseif v isa Dict{Tuple{Tuple{String,String},Tuple{String,String}}}
            param = RawParam(ks,Vector{NTuple{4,String}}(undef,0),Vector{valtype(v)}(undef,0),String[],String[],assocdata,:unknown)
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

function _fill_sources!(input,allsources,tofill)
    for i in eachindex(input)
        input[i] = coalesce(allsources[i],tofill)
    end
    return input
end

function build_raw_param(name,comps,vals,sources,csv,csvtype,grouptype)
    s::Vector{Int} = findall(!ismissing,vals)
    ls = length(s)
    _vals = Vector{nonmissingtype(eltype(vals))}(undef,ls)
    _sources = Vector{String}(undef,ls)
    _comps = similar(comps,ls)
    _csv = Vector{String}(undef,ls)
    for (i,j) ∈ pairs(s)
        _comps[i] = comps[j]
        _vals[i] = vals[j]
        _sources[i] = sources[j]
        _csv[i] = csv[j]
    end
    return RawParam(string(strip(name)),_comps,_vals,_sources,_csv,csvtype,grouptype)
end
#verbose functionality, is executed for each csv when verbose == true

function __verbose_findparams_invaliddata(filepath)
    @warn "Skipping $filepath, cannot infer correct csv type. Check line 2 of the CSV to see if it has valid information."
end

function __assoc_string(pair)
    "($(pair[1]),$(pair[3])) ⇋ ($(pair[2]), $(pair[4]))"
end

function __verbose_findparams_start(filepath,components,headerparams,parsegroups,csvtype,grouptype)
    csv_string = Symbol(csvtype)
    no_parsegroups = parsegroups == :off
    if no_parsegroups
        if csvtype ∈ (groupdata,structgroupdata)
            @info("Skipping $csv_string csv $filepath")
        else
            @info("Searching for $csv_string headers $headerparams for query $components at $filepath ...")
        end
    else
        if csvtype == groupdata
            @info("Searching for groups for components $components at $filepath ...")
        elseif csvtype == structgroupdata
            @info("Searching for intragroup interactions for components $components at $filepath ...")
        else
            @info("Skipping $csv_string csv $filepath")
        end
    end
    if grouptype != :unknown
        @info("group type: $grouptype")
    end
end


function __verbose_findparams_found(foundvalues)
    for v ∈ foundvalues
        if v.type == singledata
            io = IOBuffer()
            show_pairs(io,first.(v.component_info),v.data," => ",quote_string = false)
            vals = String(take!(io))
            kk = info_color(v.name)
            @info("""Found single component data: $kk with values:
            $vals
            """)
        elseif v.type == pairdata
            io = IOBuffer()
            show_pairs(io,first.(v.component_info,2),v.data," => ",quote_string = false)
            vals = String(take!(io))
            kk = info_color(v.name)
            @info("""Found pair component data: $kk with values:
            $vals
            """)
        elseif v.type == assocdata
            io = IOBuffer()
            show_pairs(io,__assoc_string.(v.component_info),v.data," => ",quote_string = false)
            vals = String(take!(io))
            kk = info_color(v.name)
            @info("""Found association component data: $kk with values:
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
        elseif v.type == structgroupdata
            @info("TODO: parse intragroup data for debug")
        end
    end
end

const readcsvtype_keywords  = ["like", "single", "unlike", "pair", "assoc", "association", "group", "groups","intragroup","intragroups"]

function read_csv_options(filepath)
    return _read_csv_options(getline(String(filepath), 2))
end

function _read_csv_options(line::String)
    re = r"\[.*\]"
    maybe_opts = match(re,line)

    # Searches for type from second line of CSV.
    has_csv_options = !isnothing(maybe_opts)
    if has_csv_options
        opts = chop(maybe_opts.match,head = 1,tail = 1)
        return __get_options(opts)
    else
        keywords = readcsvtype_keywords
        words = split(lowercase(strip(line, ',')), ' ')
        foundkeywords = intersect(words, keywords)
        _species = intersect(words,["species"])
        _estimator = intersect(words,["method"])
        return (csvtype = _readcsvtype(foundkeywords),grouptype = :unknown,estimator = _estimator, species = _species,sep = :comma)
    end
end

const NT_CSV_OPTIONS = (csvtype = namedtupledata,grouptype = :unknown,estimator = :no_estimator, species = ["all"],sep = :comma)
+
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
    key == "intragroup" && return structgroupdata
    key == "intragroups" && return structgroupdata
    key == "invalid" && return invaliddata
    return invaliddata
end

function __get_options(data)
    opts = eachsplit(data,',')
    opts_dict = Dict{String,String}()
    for opt in opts
        k,v = _parse_kv(opt,"=")
        opts_dict[k] = v
    end
    _csvtype = _readcsvtype(get(opts_dict,"csvtype","invalid"))
    _grouptype = Symbol(get(opts_dict,"grouptype","unknown"))
    _estimator = Symbol(get(opts_dict,"method","error"))
    _estimator = Symbol(get(opts_dict,"method","error"))
    _species = String.(split(get(opts_dict,"species","all")," "))
    _sep = Symbol(get(opts_dict,"sep","comma"))
    return (csvtype = _csvtype,grouptype = _grouptype,estimator = _estimator, species = _species, sep = _sep)
end

function valid_headerparams(csvheaders, options::ParamOptions = DefaultOptions)
    ignorelist = deepcopy(options.ignore_headers)
    push!(ignorelist,options.species_columnreference)
    push!(ignorelist,options.source_columnreference)
    push!(ignorelist,options.site_columnreference)
    map!(normalisestring,ignorelist,ignorelist)
    result = filter(csvheaders) do header
        norm_header = normalisestring(header)
        #that regex allows to ignore things like "species1" or "SMILES1"
        normalisestring(norm_header; tofilter=r"[ \-\_\d]") ∉ ignorelist
    end
end

include("database_group.jl")
include("database_writer.jl")

export getparams