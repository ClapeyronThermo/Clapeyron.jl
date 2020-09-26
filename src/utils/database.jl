include("csv_parser.jl")
import Combinatorics

function methods_return(selected_method::String = "None")
    # Returns an Array of all methods in the database directory if input is "None"
    # If input is any other string, it will return an Array of just that string
    if selected_method == "None"
        return filter(x-> isdir(joinpath(dirname(pathof(JuliaSAFT)), "../database",x)), readdir(joinpath(dirname(pathof(JuliaSAFT)), "../database")))
    else
        if selected_method in filter(x-> isdir(joinpath(dirname(pathof(JuliaSAFT)), "../database",x)), readdir(joinpath(dirname(pathof(JuliaSAFT)), "../database")))
            return [selected_method]
        else
            error("Database for selected method " * selected_method * " not found.")
        end
    end
end

function searchdatabase_like(components::Array{String, 1}, selected_method="None", user_input::AbstractString="None")
    # Returns a dictionary of components containing another dictionary
    # where the key is a method in which its database contains that parameter,
    # with the found line number as value
    methods = methods_return(selected_method)
    found_methods = Dict(Set([component]) => Dict{String, Int64}() for component in components)
    for component in components
        for method in methods
            filepath = joinpath(dirname(pathof(JuliaSAFT)), "../database", method, "data_" * method * "_like" * ".csv")
            if isfile(filepath)
                found_lines = findmatches(filepath, component, "species"; header_row=3)
                if length(found_lines) > 1
                    println("The species " * component * " is not unique in database for method " * method * ". Selecting most recent entry.")
                elseif length(found_lines) == 0
                    if length(found_lines) == 0
                        if selected_method != "None"
                            error("The species " * component * " cannot be found in database for method " * method * ".")
                        end
                    end
                end
                if length(found_lines) != 0
                    found_methods[Set([component])][method] = found_lines[end]
                end
            end
        end
    end
    return found_methods
end

function searchdatabase_unlike(components::Array{String, 1}, selected_method="None", user_input::AbstractString="None")
    # Returns a dictionary of pairs (no same pair) containing another dictionary
    # where the key is a method in which its database contains that parameter,
    # with the found line number as value
    methods = methods_return(selected_method)
    pairs = [Set(i) for i in collect(Combinatorics.combinations(components, 2))]
    found_methods = Dict{Set{String}, Dict{String, Int64}}() 
    for pair in pairs
        for method in methods
            filepath = joinpath(dirname(pathof(JuliaSAFT)), "../database", method, "data_" * method * "_unlike" * ".csv")
            if isfile(filepath)
                components = [i for i in pair]
                found_lines = findmatches_pair(filepath, components[1], components[2], "species1", "species2"; header_row=3)
                if length(found_lines) > 1
                    println("The pair (" * join([i for i in pair], ", ") * ") is not unique in database for method " * method * ". Selecting most recent entry.")
                end
                if length(found_lines) != 0
                    found_methods[pair] = Dict()
                    found_methods[pair][method] = found_lines[end]
                end
            end
        end
    end
    return found_methods
end

function searchdatabase_assoc(components::Array{String, 1}, selected_method="None", user_input::AbstractString="None")
    # Returns a dictionary of pairs (same pair allowed) containing another dictionary
    # where the key is a method in which its database contains that parameter,
    # with the found line number as value
    methods = methods_return(selected_method)
    pairs = vcat([Set([i]) for i in components], [Set(i) for i in collect(Combinatorics.combinations(components, 2))])
    found_methods = Dict{Set{String}, Dict{String, Array{Int64,1}}}() 
    for pair in pairs
        for method in methods
            filepath = joinpath(dirname(pathof(JuliaSAFT)), "../database", method, "data_" * method * "_assoc" * ".csv")
            if isfile(filepath)
                if length(pair) == 1
                    component = [i for i in pair]
                    components = [component component]
                else
                    components = [i for i in pair]
                end
                found_lines = findmatches_pair(filepath, components[1], components[2], "species1", "species2"; header_row=3)
                if length(found_lines) > 0
                    found_methods[pair] = Dict()
                    found_methods[pair][method] = found_lines
                end
            end
        end
    end
    return found_methods
end

function retrieveparams_like(components::Array{String, 1}, selected_method, user_input::AbstractString="None")
    # Returns a dictionary of all columns in the like database for a selected method
    # for each component
    filepath = joinpath(dirname(pathof(JuliaSAFT)), "../database", selected_method, "data_" * selected_method * "_like" * ".csv")
    header = parseline(filepath, 3)
    found_method = searchdatabase_like(components, selected_method)
    found_params = Dict{Set{String}, Dict{String, Any}}() 
    for component in keys(found_method)
        found_params[component] = Dict(zip(header, parseline(filepath, found_method[component][selected_method])))
    end
    return found_params
end

function retrieveparams_unlike(components::Array{String, 1}, selected_method, user_input::AbstractString="None")
    # Returns a dictionary of all columns in the unlike database for a selected method
    # for each pair
    filepath = joinpath(dirname(pathof(JuliaSAFT)), "../database", selected_method, "data_" * selected_method * "_unlike" * ".csv")
    header = parseline(filepath, 3)
    found_method = searchdatabase_unlike(components, selected_method)
    found_params = Dict{Set{String}, Dict{String, Any}}() 
    for pair in keys(found_method)
        found_params[pair] = Dict(zip(header, parseline(filepath, found_method[pair][selected_method])))
    end
    return found_params
end

function retrieveparams_assoc(components::Array{String, 1}, selected_method, user_input::AbstractString="None")
    # Returns a dictionary of all columns in the assoc database for a selected method
    # for each pair (same pair allowed)
    filepath = joinpath(dirname(pathof(JuliaSAFT)), "../database", selected_method, "data_" * selected_method * "_assoc" * ".csv")
    header = parseline(filepath, 3)
    found_method = searchdatabase_assoc(components, selected_method)
    found_params = Dict{Set{String}, Dict{String, Dict{Set{String}, Dict{String, Any}}}}() 
    if !isempty(found_method)
        for pair in keys(found_method)
            found_params[pair] = Dict("assoc" => Dict())
            for line_number in found_method[pair][selected_method]
                retrieved = Dict(zip(header, parseline(filepath, line_number)))
                assoc_pair = Set([retrieved["site1"], retrieved["site2"]])
                found_params[pair]["assoc"][assoc_pair] = retrieved 
            end
        end
    end
    return found_params
end

function retrieveparams(components::Array{String, 1}, selected_method, user_input::AbstractString="None")
    # Combines all retrieved parameters into a single dictionary
    # The keys will be sets of single and paired components
    params_like = retrieveparams_like(components, selected_method)
    params_unlike = retrieveparams_unlike(components, selected_method)
    params_assoc = retrieveparams_assoc(components, selected_method)
    pairs = union(Set(keys(params_like)), Set(keys(params_unlike)), Set(keys(params_assoc)))
    params = Dict{Set{String}, Any}()
    for pair in pairs
        params[pair] = Dict()
        if haskey(params_like, pair)
            merge!(params[pair], params_like[pair])
        end
        if haskey(params_unlike, pair)
            merge!(params[pair], params_unlike[pair])
        end
        if haskey(params_assoc, pair)
            merge!(params[pair], params_assoc[pair])
        end
    end
    return params
end

function filterparams(raw_params::Dict{Set{String}}, pure_params::T; pair_params::T=Array{String,1}([]), assoc_params::T=Array{String,1}([])) where T<:Array{String,1}
    # Filters the raw parametrs from retrieveparameters into selected headers
    # Returns dictionaries where the keys are the header of teh selected columns
    # One dictionary for each like, unlike, and assoc
    components = filter(x -> length(x)==1, keys(raw_params))
    pairs = filter(x -> length(x)==2, keys(raw_params))
    pure_params_dict = Dict{String, Dict{Set{String}, Float64}}()
    pair_params_dict = Dict{String, Dict{Set{String}, Float64}}()
    assoc_params_dict = Dict{String, Dict{Set{String}, Dict{Set{String},Float64}}}()
    for pure_param in pure_params
        pure_params_dict[pure_param] = Dict()
        for component in components 
            param_value = raw_params[component][pure_param]
            if !ismissing(param_value)
                push!(pure_params_dict[pure_param], component => param_value)
            end
        end
    end
    for pair_param in pair_params
        pair_params_dict[pair_param] = Dict()
        for pair in pairs
            if haskey(raw_params[pair], pair_param)
                param_value = raw_params[pair][pair_param]
                if !ismissing(param_value)
                    push!(pair_params_dict[pair_param], pair => param_value)
                end
            end
        end
    end
    for assoc_param in assoc_params
        assoc_params_dict[assoc_param] = Dict()
        for component in union(components, pairs)
            if haskey(raw_params[component], "assoc")
                for assoc_pair in keys(raw_params[component]["assoc"])
                    param_value = raw_params[component]["assoc"][assoc_pair][assoc_param]
                    if !ismissing(param_value)
                        if !haskey(assoc_params_dict[assoc_param], component)
                            assoc_params_dict[assoc_param][assoc_pair] = Dict()
                        end
                        push!(assoc_params_dict[assoc_param][assoc_pair], component => param_value)
                    end
                end
            end
        end
    end
    #= for assoc_param in assoc_params =#
    #=     assoc_params_dict[assoc_param] = Dict() =#
    #=     for component in union(components, pairs) =#
    #=         if haskey(raw_params[component], "assoc") =#
    #=             for assoc_pair in keys(raw_params[component]["assoc"]) =#
    #=                 param_value = raw_params[component]["assoc"][assoc_pair][assoc_param] =#
    #=                 if !ismissing(param_value) =#
    #=                     if !haskey(assoc_params_dict[assoc_param], component) =#
    #=                         assoc_params_dict[assoc_param][component] = Dict() =#
    #=                     end =#
    #=                     push!(assoc_params_dict[assoc_param][component], assoc_pair => param_value) =#
    #=                 end =#
    #=             end =#
    #=         end =#
    #=     end =#
    #= end =#
    return pure_params_dict, pair_params_dict, assoc_params_dict
end
    
