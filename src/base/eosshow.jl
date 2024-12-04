"""
custom_show(::Type{T})::Bool

Trait used to determine if an EoSModel can use the default custom `show` methods defined for EoSModels in the package.

"""
custom_show(model::EoSModel) = _custom_show(model)
custom_show(model) = false
function _custom_show(Base.@nospecialize(model))
    hasfield(typeof(model),:components)
end 

#function used to customize the first line to your liking
show_info(io,model) = nothing

function show_params(io,model)
    hasfield(typeof(model),:params) || return nothing
    iszero(fieldcount(typeof(model.params))) && return nothing
    println(io)
    paramnames = fieldnames(typeof(model.params))
    len_params = length(paramnames)
    !iszero(len_params) && print(io,"Contains parameters: ")
    show_pairs(io,paramnames,pair_separator = ", ",quote_string = false)
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

function eosshow(io::IO, mime::MIME"text/plain", Base.@nospecialize(model::EoSModel))   
    print(io, typeof(model))
    if hasfield(typeof(model),:components)
        length(model) == 1 && println(io, " with 1 component:")
        length(model) > 1 && println(io, " with ", length(model), " components:")
        if has_groups(model)
            groups = model.groups
            show_groups(io,groups)
            println(io)
            print(io,"Group Type: ",groups.grouptype)
        else
            show_pairs(io,model.components)
        end
    else
        print(io,"()")
    end
    show_info(io,model)
    show_params(io,model)
    show_reference_state(io,model)
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

#overload of Base.show here
function Base.show(io::IO,mime::MIME"text/plain",model::EoSModel)
    if custom_show(model)
        eosshow(io,mime,model)
    else
        show_default(io,mime,model)
    end
end

function Base.show(io::IO,model::EoSModel)
    if custom_show(model)
        eosshow(io,model)
    else
        show_default(io,model)
    end
end

function show_reference_state(io::IO,model::EoSModel)
    return show_reference_state(io,reference_state(model),model)
end

show_reference_state(io::IO,ref::Nothing,model::EoSModel) = nothing

function show_reference_state(io::IO,ref,model::EoSModel)
    type = ref.std_type
    if type != :no_set
        println(io)
        print(io,"Reference state: ")
        print(io,type)
    end
end

export eosshow