PARAM_LOCATION = ""
import UUIDs
const PKG_UUID = parse(UUIDs.UUID,"7c7805af-46cc-48c9-995b-ed0ed2dc909a")

#not constant,as is modified at runtime
#performance is not priority

function random_csv_name(table_type)
    #return repr(table_type) * repr(rand(UInt))
    return "$(Symbol(table_type))_$(repr(rand(UInt))).csv"
end
"""
    ParamTable(type::Symbol,table;
    location = nothing,
    name = nothing,
    grouptype = :unknown,
    options = ParamOptions())

Creates a clapeyron CSV file and returns the location of that file. the type determines the table type:
- `:single` creates a table with single parameters
- `:pair` creates a table with pair parameters
- `:assoc` creates a table with association parameters
- `:group` creates a table with association parameters
By default, the name is generated randomly, and the table is stored as a temporary scratch space (provided by Scratch.jl).
You can clean said scratch space by using `Clapeyron.cleartemp!()`.

## Examples:
```julia-repl
julia> data = (species = ["water"],Mw = [18.03]) #it could be a Dict, a named tuple, or any Tables.jl compatible table
(species = ["water"], Mw = [18.9])
julia> file = ParamTable(:single,data ,name="water_new_mw")
"C:\\Users\\user\\.julia\\scratchspaces\\7c7805af-46cc-48c9-995b-ed0ed2dc909a\\ParamTables\\singledata_water_new_mw.csv"
julia> model = PCSAFT(["water","methanol"],userlocations = [file])
PCSAFT{BasicIdeal} with 2 components:
 "water"
 "methanol"
Contains parameters: Mw, segment, sigma,
epsilon, epsilon_assoc, bondvol
julia> model.params.Mw
SingleParam{Float64}("Mw") with 2 components:
 "water" => 18.9
 "methanol" => 32.042
```
"""
function ParamTable(type::Symbol,data;
    location::Union{String,Nothing} = nothing,
    name::Union{String,Nothing} = nothing,
    grouptype::Symbol = :unknown,
    options::ParamOptions = DefaultOptions)
    if location === nothing
        location = generate_location!()
    end
    table_type = _readcsvtype(String(type))
    if name === nothing
        name = repr(rand(UInt))
    end
    csvname = "$(Symbol(table_type))_$(name).csv"
    headers = String.(colnames(data))
    normalised_headers = normalisestring.(headers)
    _,_,_ = col_indices(table_type,normalised_headers,options) #basically to check the schema
    file = joinpath(location,csvname)
    io = open(file,"w")
    pretext =
    """Clapeyron Database File
    $(name) Parameters [csvtype = $type, grouptype = $grouptype]
    """
    write(io,pretext)
    CSV.write(io, data,append = true,header = true)
    close(io)
    return file
end

function generate_location!()
    if PARAM_LOCATION === ""
        global PARAM_LOCATION = @get_scratch!("ParamTables")
    end
    return PARAM_LOCATION
end
"""
    cleartemp!()
Deletes all files in the temporary Clapeyron scratch space, used to store the csvs created by `ParamTable`.
"""
function cleartemp!()
    Scratch.delete_scratch!(PKG_UUID,"ParamTables")
end

colnames(x) = Tables.columnnames(x)

export ParamTable