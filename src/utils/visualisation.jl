

function eosshow(io::IO, ::MIME"text/plain", Base.@nospecialize(model::EoSModel))
    print(io, typeof(model))
    length(model) == 1 && println(io, " with 1 component:")
    length(model) > 1 && println(io, " with ", length(model), " components:")
    for i in 1:length(model)
        print(io, " \"", model.components[i], "\"")
        println(io)
    end
    if hasfield(typeof(model),:params)
        paramnames = fieldnames(typeof(model.params))
        len_params = length(paramnames)
        !iszero(len_params) && print(io,"Contains parameters: ")
        firstloop = true
        for fieldname in paramnames
            firstloop == false && print(io, ", ")
            print(io, fieldname)
            firstloop = false
        end
    end
end

function eosshow(io::IO, Base.@nospecialize(model::EoSModel))
    print(io, typeof(model))
    firstloop = true
    print(io, "(")
    for i in 1:length(model.components)
        firstloop == false && print(io, ", ")
        print(io, "\"", model.components[i], "\"")
        firstloop = false
    end
    print(io, ")")
end

function gc_eosshow(io::IO, ::MIME"text/plain", Base.@nospecialize(model::EoSModel))
    print(io, typeof(model))
    length(model) == 1 && println(io, " with 1 component:")
    length(model) > 1 && println(io, " with ", length(model), " components:")
    for i in 1:length(model)
        print(io, " \"", model.components[i], "\": ")
        firstloop = true
        for k in 1:length(model.groups.groups[i])
            firstloop == false && print(io, ", ")
            print(io, "\"", model.groups.groups[i][k], "\" => ", model.groups.n_groups[i][k])
            firstloop = false
        end
        println(io)
    end
    print(io, "Contains parameters: ")
    firstloop = true
    for fieldname in fieldnames(typeof(model.params))
        firstloop == false && print(io, ", ")
        print(io, fieldname)
        firstloop = false
    end
end


export eosshow