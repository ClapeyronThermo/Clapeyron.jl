
function show_params(io,model)
    if hasfield(typeof(model),:params)
        println(io)
        paramnames = fieldnames(typeof(model.params))
        len_params = length(paramnames)
        !iszero(len_params) && print(io,"Contains parameters: ")
        show_pairs(io,paramnames,pair_separator = ", ",quote_string = false)
    end
end

function may_show_references(io::IO,model)
    if get(ENV,"CLAPEYRON_SHOW_REFERENCES","FALSE") == "TRUE"
        show_references(io,model)     
    end
end

function show_references(io::IO,model)
    citations = cite(model)
    iszero(length(citations)) && return nothing #do not do anything if there isnt any citations
    println(io)
    print(io,"References: ") 
    for (i,doi) in enumerate(cite(model))
        i != 1 && print(io,", ")
        print(io,doi)   
    end
end


function eosshow(io::IO, ::MIME"text/plain", Base.@nospecialize(model::EoSModel))
    print(io, typeof(model))
    if hasfield(typeof(model),:components)
        length(model) == 1 && println(io, " with 1 component:")
        length(model) > 1 && println(io, " with ", length(model), " components:")
        show_pairs(io,model.components)
    else
        print(io,"()")
    end

    show_params(io,model)
    may_show_references(io,model)
end

function eosshow(io::IO, Base.@nospecialize(model::EoSModel))
    print(io, typeof(model))
    print(io, "(")
    show_pairs(io,model.components,pair_separator = ", ")
    print(io, ")")
end

#utilities for showing groups
function __show_group_i(io,val,missingvalue = "")
    keys,vals = val
    #@show val
    if length(vals) == 0 && missingvalue != ""
        print(io,missingvalue)
    else
        show_pairs(io,keys,vals," => ",pair_separator = ", ")
    end
end

show_groups(io,gc) = show_pairs(io,gc.components,zip(gc.groups,gc.n_groups),": ",__show_group_i)

function gc_eosshow(io::IO, ::MIME"text/plain", Base.@nospecialize(model::EoSModel))
    print(io, typeof(model))
    length(model) == 1 && println(io, " with 1 component:")
    length(model) > 1 && println(io, " with ", length(model), " components:")
    groups = model.groups
    show_groups(io,groups)
    println(io)
    print(io,"Group Type: ",groups.grouptype)
    show_params(io,model)
    may_show_references(io,model)
end

function gc_eosshow(io::IO, Base.@nospecialize(model::EoSModel))
    return eosshow(io,model)
end

export eosshow