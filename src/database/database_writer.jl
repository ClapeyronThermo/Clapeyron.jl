function get_header(path)
    io = IOBuffer()
    _sep = Any[1]
    file = open(getpath(path),"r")
    for (k,i) in enumerate(eachline(file))
        
        if k == 2
            io2 = IOBuffer()
            print(io2,i)
            str = String(take!(io2))
            _sep[] = _read_csv_options(str)
        end
        if k == 3
            break
        end
        println(io,i)
    end
    header = String(take!(io))
    close(file)
    _delims = (comma = ',',space = ' ')
    sep = _sep[][:sep]
    if sep isa Symbol
        _delim = get(_delims,sep,string(sep))
    else
        _delim = sep
    end
    return header,_delim,_sep[]
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

function make_header(_database,name = nothing,grouptype = :unknown)
    
    io = IOBuffer()
    println(io,"Clapeyron Database File")
    if name !== nothing
        print(io,name)
    end

    #create a csv options type, based on the column names stored.
    headers = normalisestring.(String.(Tables.columnnames(_database)))
    if all(in(headers),("species","groups"))
        csvtype = :group
    elseif all(in(headers),("species1","species2","site1","site2"))
        csvtype = :assoc
    elseif all(in(headers),("species1","species2"))
        csvtype = :pair
    elseif all(in(headers),("species"))
        csvtype = :like
    else
        csvtype = :invalid
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
    sep = ','
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





