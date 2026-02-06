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

"""
    eosshow(io::IO, model::EoSModel)
    eosshow(io::IO, ::MIME"text/plain", model::EoSModel)

Custom pretty-printer for `EoSModel` instances.

This is the backend used by `Base.show` for models that opt into the custom
display. The text/plain variant prints components, parameters, reference state,
and (optionally) citations when enabled via `ENV["CLAPEYRON_SHOW_REFERENCES"]`.
"""
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
            show_pairs(io,component_list(model))
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
    show_pairs(io,component_list(model),pair_separator = ", ")
    print(io, ")")
end

#utilities for showing groups
function __show_group_ij(io,v)
    if isinteger(v)
        print(io,Int(v))
    else
        print(io,v)
    end
end

function __show_group_i(io,val,missingvalue = "")
    keys,vals = val
    #@show val
    if length(vals) == 0 && missingvalue != ""
        print(io,missingvalue)
    else
        show_pairs(io,keys,vals," => ",pair_separator = ", ",__show_group_ij)
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

function show_reference_state(io::IO,model::EoSModel;space = false)
    return show_reference_state(io,reference_state(model),model,space)
end

show_reference_state(io::IO,ref::Nothing,model::EoSModel,space) = nothing

function show_reference_state(io::IO,ref,model::EoSModel,space)
    type = ref.std_type
    if type != :no_set
        println(io)
        space && print(io," ")
        print(io,"Reference state: ")
        print(io,type)
    end
end

"""
    eos_repr(io::IO,model,newlines = true)
    eos_repr(model;newlines = true)::String

Given a `model::EoSModel`, returns a copy-pastable text representation of that model, designed to be parseable julia code.
By default, `eos_repr` inserts some newlines for easier human reading of the output, that can be disabled via the `newlines` keyword argument.

## Examples:

```julia-repr
julia> model = PCSAFT("methane")
PCSAFT{BasicIdeal, Float64} with 1 component:
 "methane"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> push!(model.references,"test")
3-element Vector{String}:
 "10.1021/ie0003887"
 "10.1021/ie010954d"
 "test"

julia> s = eos_repr(model);

julia> model2 = eval(Meta.parse(s)) #exactly the same model
PCSAFT{BasicIdeal, Float64} with 1 component:
 "methane"
Contains parameters: Mw, segment, sigma, epsilon, epsilon_assoc, bondvol

julia> model2.references
3-element Vector{String}:
 "10.1021/ie0003887"
 "10.1021/ie010954d"
 "test"
```
"""
function eos_repr(io::IO,model;inside = false,newlines = true)
    M = typeof(model)
    n = fieldnames(M)
    newlines && print(io,"\n")
    print(io,M)
    print(io,"(")
    if !inside && newlines
        print(io,"\n")
    end
    k = 0
    nf = length(n)
    for i in n
        k += 1
        f = getfield(model,i)
        F = typeof(f)
        if F.name.module == Clapeyron || F isa EoSModel
            eos_repr(io,f,inside = true)
        elseif f isa PackedVofV
            print(io,"Clapeyron.PackedVofV(")
            show(io,f.p)
            print(io,", ")
            show(io,f.v)
            print(io,")")
        else
            show(io,f)
        end
        k != nf && print(io,", ")
    end
    print(io,")")


    return nothing
end

function eos_repr(model;newlines = false)
    io = IOBuffer()
    eos_repr(io,model;newlines)
    return String(take!(io))
end

export eosshow, eos_repr