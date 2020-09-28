include("csv_parser.jl")

function methods_return(selected_method::String = "None")
    # Returns an Array of all methods in the database directory if input is "None"
    # If input is any other string, it will return an Array of just that string
    if selected_method == "None"
        return filter(x-> isdir(joinpath(dirname(pathof(JuliaSAFT)), "../database", x)), readdir(joinpath(dirname(pathof(JuliaSAFT)), "../database")))
    else
        if selected_method in filter(x-> isdir(joinpath(dirname(pathof(JuliaSAFT)), "../database",x)), readdir(joinpath(dirname(pathof(JuliaSAFT)), "../database")))
            return [selected_method]
        else
            error("Database for selected method $selected_method not found.")
        end
    end
end

function createfilepath(selected_method, database_type; variant="None")
    return joinpath(dirname(pathof(JuliaSAFT)), "../database", selected_method, variant=="None" ? "data_$(selected_method)_$(database_type).csv" : "$variant/data_$(selected_method)_$(variant)_$(database_type).csv")
end

function customdatabase_check(selected_method, database_type, customdatabase_filepath; variant="None")
    if !isfile(customdatabase_filepath)
        error("Custom database $customdatabase_filepath not found.")
    elseif selected_method == "None"
        error("Method has to be specified for custom databases.")
    end
    database_filepath = createfilepath(selected_method, database_type; variant=variant)
    database_header = parseline(database_filepath, 3)
    customdatabase_header = parseline(customdatabase_filepath, 3)
    if Set(database_header) != Set(customdatabase_header)
        error("Custom database $customdatabase_filepath header not consistent with format in $selected_method (variant: $variant) for type $database_type.")
    end
end


function searchdatabase_like(components::Array{String, 1}, selected_method="None"; customdatabase_filepath="None", variant="None")
    # Returns a dictionary of components containing another dictionary
    # where the key is a method in which its database contains that parameter,
    # with the found line number as value
    if customdatabase_filepath != "None"
        customdatabase_check(selected_method, "like", customdatabase_filepath, variant=variant)
    end
    methods = methods_return(selected_method)
    found_methods = Dict(Set([component]) => Dict{String, Int64}() for component in components)
    for component in components
        for method in methods
            filepath = customdatabase_filepath == "None" ? createfilepath(selected_method, "like"; variant=variant) : customdatabase_filepath
            if isfile(filepath)
                found_lines = findmatches(filepath, component, "species"; header_row=3)
                if length(found_lines) > 1
                    println("The species $component is not unique in database for method $method. Selecting most recent entry.")
                elseif length(found_lines) == 0
                    if length(found_lines) == 0
                        if selected_method != "None"
                            error("The species $component cannot be found in database for method $method.")
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

function searchdatabase_unlike(components::Array{String, 1}, selected_method="None"; customdatabase_filepath="None", variant="None")
    # Returns a dictionary of pairs (no same pair) containing another dictionary
    # where the key is a method in which its database contains that parameter,
    # with the found line number as value
    if customdatabase_filepath != "None"
        customdatabase_check(selected_method, "unlike", customdatabase_filepath, variant=variant)
    end
    methods = methods_return(selected_method)
    pairs = [Set(i) for i in collect(combinations(components, 2))]
    found_methods = Dict{Set{String}, Dict{String, Int64}}()
    for pair in pairs
        for method in methods
            filepath = customdatabase_filepath == "None" ? createfilepath(selected_method, "unlike"; variant=variant) : customdatabase_filepath
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

function searchdatabase_assoc(components::Array{String, 1}, selected_method="None"; customdatabase_filepath="None", variant="None")
    # Returns a dictionary of pairs (same pair allowed) containing another dictionary
    # where the key is a method in which its database contains that parameter,
    # with the found line number as value
    if customdatabase_filepath != "None"
        customdatabase_check(selected_method, "assoc", customdatabase_filepath, variant=variant)
    end
    methods = methods_return(selected_method)
    pairs = vcat([Set([i]) for i in components], [Set(i) for i in collect(Combinatorics.combinations(components, 2))])
    found_methods = Dict{Set{String}, Dict{String, Array{Int64,1}}}()
    for pair in pairs
        for method in methods
            filepath = customdatabase_filepath == "None" ? createfilepath(selected_method, "assoc"; variant=variant) : customdatabase_filepath
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

function retrieveparams_like(components::Array{String, 1}, selected_method; customdatabase_filepath="None", variant="None", redirect="None")
    # Returns a dictionary of all columns in the like database for a selected method
    # for each component
    filepath = customdatabase_filepath == "None" ? createfilepath(selected_method, "like"; variant=variant) : customdatabase_filepath
    if redirect != "None"
        selected_method = redirect
    end
    header = parseline(filepath, 3)
    found_method = searchdatabase_like(components, selected_method; customdatabase_filepath=customdatabase_filepath, variant=variant)
    found_params = Dict{Set{String}, Dict{String, Any}}()
    for component in keys(found_method)
        found_params[component] = Dict(zip(header, parseline(filepath, found_method[component][selected_method])))
    end
    return found_params
end

function retrieveparams_unlike(components::Array{String, 1}, selected_method; customdatabase_filepath="None", variant="None", redirect="None")
    # Returns a dictionary of all columns in the unlike database for a selected method
    # for each pair
    filepath = customdatabase_filepath == "None" ? createfilepath(selected_method, "unlike"; variant=variant) : customdatabase_filepath
    if redirect != "None"
        selected_method = redirect
    end
    header = parseline(filepath, 3)
    found_method = searchdatabase_unlike(components, selected_method; customdatabase_filepath=customdatabase_filepath, variant=variant)
    found_params = Dict{Set{String}, Dict{String, Any}}()
    for pair in keys(found_method)
        found_params[pair] = Dict(zip(header, parseline(filepath, found_method[pair][selected_method])))
    end
    return found_params
end

function retrieveparams_assoc(components::Array{String, 1}, selected_method; customdatabase_filepath="None", variant="None", redirect="None")
    # Returns a dictionary of all columns in the assoc database for a selected method
    # for each pair (same pair allowed)
    filepath = customdatabase_filepath == "None" ? createfilepath(selected_method, "assoc"; variant=variant) : customdatabase_filepath
    if redirect != "None"
        selected_method = redirect
    end
    header = parseline(filepath, 3)
    found_method = searchdatabase_assoc(components, selected_method; customdatabase_filepath=customdatabase_filepath, variant=variant)
    pairs = keys(found_method)
    found_params = Dict{Set{String}, Dict{Set{Tuple{String, String}}, Dict{String, Any}}}()
    if !isempty(found_method)
        for pair in pairs
            found_params[pair] = Dict()
            for line_number in found_method[pair][selected_method]
                retrieved = Dict(zip(header, parseline(filepath, line_number)))
                assoc_pair = Set([(retrieved["species1"], retrieved["site1"]), (retrieved["species2"], retrieved["site2"])])
                found_params[pair][assoc_pair] = retrieved
            end
        end
        # Check if reverse pair exists; if not, make create reverse pair equal to original pair
        for pair in pairs
            for assoc_pair in keys(found_params[Set(pair)])
                if length(assoc_pair) == 2
                    pair_tuples = collect(assoc_pair)
                    reverse_assoc_pair = Set([(pair_tuples[1][1], pair_tuples[2][2]), (pair_tuples[2][1], pair_tuples[1][2])])
                    if !haskey(found_params[Set(pair)], reverse_assoc_pair)
                        found_params[Set(pair)][reverse_assoc_pair] =
                            found_params[Set(pair)][assoc_pair]
                    end
                end
            end
        end
    end
    return found_params
end

function retrieveparams(components::Array{String, 1}, selected_method;
        customdatabase_like = "None", variant_like = "None", redirect_like = "None",
        customdatabase_unlike = "None", variant_unlike = "None", redirect_unlike = "None",
        customdatabase_assoc = "None", variant_assoc = "None", redirect_assoc = "None")
    params_like = retrieveparams_like(components, selected_method; customdatabase_filepath = customdatabase_like, variant = variant_like, redirect = redirect_like)
    params_unlike = retrieveparams_unlike(components, selected_method; customdatabase_filepath = customdatabase_unlike, variant = variant_unlike, redirect = redirect_unlike)
    params_assoc = retrieveparams_assoc(components, selected_method; customdatabase_filepath = customdatabase_assoc, variant = variant_assoc, redirect = redirect_assoc)
    return [params_like, params_unlike, params_assoc]
end

function filterparams(raw_params, like_params::T; unlike_params::T=Array{String,1}([]), assoc_params::T=Array{String,1}([])) where T<:Array{String,1}
    # Filters the raw parametrs from retrieveparameters into selected headers
    # Returns dictionaries where the keys are the header of the selected columns
    # One dictionary for each like, unlike, and assoc
    raw_params_like = raw_params[1]
    raw_params_unlike = raw_params[2]
    raw_params_assoc = raw_params[3]
    components = keys(raw_params_like)
    pairs = keys(raw_params_unlike)
    pair_tuples = keys(raw_params_assoc)
    like_params_dict = Dict{String, Dict{Set{String}, Float64}}()
    unlike_params_dict = Dict{String, Dict{Set{String}, Float64}}()
    assoc_params_dict = Dict{String, Dict{Set{String}, Dict{Set{Tuple{String, String}}, Float64}}}()
    for pure_param in like_params
        like_params_dict[pure_param] = Dict()
        for component in components
            param_value = raw_params_like[component][pure_param]
            if !ismissing(param_value)
                push!(like_params_dict[pure_param], component => param_value)
            end
        end
    end
    for pair_param in unlike_params
        unlike_params_dict[pair_param] = Dict()
        for pair in pairs
            if haskey(raw_params_unlike[pair], pair_param)
                param_value = raw_params_unlike[pair][pair_param]
                if !ismissing(param_value)
                    push!(unlike_params_dict[pair_param], pair => param_value)
                end
            end
        end
    end
    for assoc_param in assoc_params
        assoc_params_dict[assoc_param] = Dict()
        for pair_tuple in pair_tuples
            for assoc_pair in keys(raw_params_assoc[pair_tuple])
                param_value = raw_params_assoc[pair_tuple][assoc_pair][assoc_param]
                if !ismissing(param_value)
                    if !haskey(assoc_params_dict[assoc_param], pair_tuple)
                        assoc_params_dict[assoc_param][pair_tuple] = Dict()
                    end
                    push!(assoc_params_dict[assoc_param][pair_tuple], assoc_pair => param_value)
                end
            end
        end
    end
    return like_params_dict, unlike_params_dict, assoc_params_dict
end
