
"""
    get_header(path)

Given a path to a Clapeyron CSV file, `get_header` will return the text of the first 3 columns, the string used as a separator, and a named tuple with all available Clapeyron CSV options.

"""
function get_header(path)
    io = IOBuffer()
    _sep = Any[1]
    file = open(getpath(path),"r")
    local str_options::String
    for (k,i) in enumerate(eachline(file)) 
        if k == 2
            io2 = IOBuffer()
            print(io2,i)
            str_options = String(take!(io2))
        end
        k == 3 && break
        println(io,i)
    end
    _sep = _read_csv_options(str_options)
    header = String(take!(io))
    close(file)
    _delims = (comma = ",",space = " ")
    sep = _sep[:sep]
    if sep isa Symbol
        _delim = get(_delims,sep,string(sep))
    else
        _delim = String(sep)
    end
    return header,_delim,_sep
end

#overwrites csv in a path
function write_csv!(path,_database;relativetodatabase = true)
    database = deepcopy(normalise_db(_database))
    header,delim,_ = get_header(path)
    io = IOBuffer()
    print(io,header)
    CSV.write(io,database,delim = delim,writeheader = true,append = true)
    normpath = getpath(path)
    f = open(normpath,"w")
    write(f,String(take!(io)))
    close(f)
    return nothing
end

"""
    make_header(_database,name = nothing,grouptype = :unknown)

Given a database, `make_header` will guess a header based on a set of rules:
- If the database only has a `species` column, then the database is considered to store single-component data
- If the database has `species1` and `species2` column names, then the database is considered to store pair-component data
- If the database has `species1`, `species2`, `site1` and `site2` then the database is considered to store association data
- If the database has `species` and `groups` then the database is considered to store group data
- If the database does not match any of the above, it is considered to store "invalid" data

It returns a string for the first 3 lines of the csv, a string for the separator, and the named tuple of Clapeyron CSV options. It has the same return type as `Clapeyron.get_header(path)`
"""
function make_header(_database,name = nothing,grouptype = :unknown)
    io = IOBuffer()
    println(io,"Clapeyron Database File")
    if name !== nothing
        print(io,name)
    end

    #create a csv options type, based on the column names stored.
    headers = Set(normalisestring.(String.(Tables.columnnames(_database))))
    if "species" in headers && "groups" in headers
        csvtype = :group
    elseif "species1" in headers && "species2" in headers && "site1" in headers && "site2" in headers
        csvtype = :assoc
    elseif "species1" in headers && "species1" in headers
        csvtype = :pair
    elseif "species" in headers

        csvtype = :like
    else
        csvtype = :invalid
    end

    use_comma = true
    if csvtype in (:group,:like)
        sp = Tables.getcolumn(_database,:species)
        for spi in sp 
            if ',' in spi 
                use_comma = false
                break
            end
        end
        sep = ifelse(use_comma,",",";")
    elseif csvtype in (:pair,:assoc)
        sp1 = Tables.getcolumn(_database,:species1)
        sp2 = Tables.getcolumn(_database,:species2)
        for (sp1i,sp2i) in zip(sp1,sp2)
            if ',' in sp1i ||  ',' in sp2i
                use_comma = false
                break
            end
        end
        sep = ifelse(use_comma,",",";")
    else
    end
    
    print(io,"[csvtype = $csvtype")
    if grouptype == :unknown
        opts = (;csvtype)
        print(io,", grouptype = $grouptype")
    else
        opts = (;csvtype,grouptype)
    end
    print(io,"]")
    header = String(take!(io))
    return header,sep,opts
end

#creates a new csv, be sure that the path actually exists.
function write_csv(path,database;relativetodatabase = true,grouptype = :unknown, name = nothing)
    _database = normalise_db(database)
    database = deepcopy(_database)
    header,delim,_ = make_header(_database,name,grouptype)
    io = IOBuffer()
    print(io,header)
    CSV.write(io,database,delim = delim,writeheader = true,append = true)
    normpath = getpath(path;relativetodatabase)
    f = open(normpath,"w")
    write(f,String(take!(io)))
    close(f)
    return nothing
end

function normalise_db(db)
    Tables.istable(typeof(db)) && return db
    throw(error("cannot normalise database into a Tables.jl compatible type."))
end

function normalise_db(db::AbstractDict{K,Any}) where {K<:AbstractString}
    return OrderedCollections.OrderedDict{Symbol,AbstractVector}(Symbol(k) => db[k] for k in keys(db))
end





