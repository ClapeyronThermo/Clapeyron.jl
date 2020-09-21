include("CSVParser.jl")
import Combinatorics

function methods_return(selected_method)
    if selected_method == "None"
        return filter(x-> isdir(joinpath("database",x)), readdir("database"))
    else
        if selected_method in filter(x-> isdir(joinpath("database",x)), readdir("database"))
            return [selected_method]
        else
            error("Database for selected method " * selected_method * " not found.")
        end
    end
end

function search_database_like(components::Array{String, 1}, selected_method="None", user_input::AbstractString="None")
    methods = methods_return(selected_method)
    found_methods = Dict(Set([component]) => Dict{String, Int64}() for component in components)
    for component in components
        for method in methods
            filepath = joinpath("database", method, "data_" * method * "_like" * ".csv")
            if isfile(filepath)
                found_lines = find_matches(filepath, component, "species"; header_row=3)
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

function search_database_unlike(components::Array{String, 1}, selected_method="None", user_input::AbstractString="None")
    methods = methods_return(selected_method)
    pairs = [Set(i) for i in collect(Combinatorics.combinations(components, 2))]
    found_methods = Dict{Set{String}, Dict{String, Int64}}() 
    for pair in pairs
        for method in methods
            filepath = joinpath("database", method, "data_" * method * "_unlike" * ".csv")
            if isfile(filepath)
                components = [i for i in pair]
                found_lines = find_matches_pair(filepath, components[1], components[2], "species1", "species2"; header_row=3)
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

function search_database_assoc(components::Array{String, 1}, selected_method="None", user_input::AbstractString="None")
    methods = methods_return(selected_method)
    pairs = vcat([Set([i]) for i in components], [Set(i) for i in collect(Combinatorics.combinations(components, 2))])
    found_methods = Dict{Set{String}, Dict{String, Array{Int64,1}}}() 
    for pair in pairs
        for method in methods
            filepath = joinpath("database", method, "data_" * method * "_assoc" * ".csv")
            if isfile(filepath)
                if length(pair) == 1
                    component = [i for i in pair]
                    components = [component component]
                else
                    components = [i for i in pair]
                end
                found_lines = find_matches_pair(filepath, components[1], components[2], "species1", "species2"; header_row=3)
                if length(found_lines) > 0
                    found_methods[pair] = Dict()
                    found_methods[pair][method] = found_lines
                end
            end
        end
    end
    return found_methods
end
    
function retrieve_parameters_assoc(components::Array{String, 1}, selected_method, user_input::AbstractString="None")
    filepath = joinpath("database", selected_method, "data_" * selected_method * "_assoc" * ".csv")
    header = read_line(filepath, 3)
    found_method = search_database_assoc(components, selected_method)
    found_parameters = Dict{Set{String}, Dict{String, Dict{Set{String}, Dict{String, Any}}}}() 
    #= found_parameters = Dict{Set{String}, Dict{String, Any}}() =# 
    if !isempty(found_method)
        for pair in keys(found_method)
            found_parameters[pair] = Dict("Assoc" => Dict())
            for line_number in found_method[pair][selected_method]
                retrieved = Dict(zip(header, tryparse_fallback(read_line(filepath, line_number))))
                assoc_pair = Set([retrieved["site1"], retrieved["site2"]])
                found_parameters[pair]["Assoc"][assoc_pair] = Dict()
                found_parameters[pair]["Assoc"][assoc_pair] = retrieved 
            end
        end
    end
    return found_parameters
end

function retrieve_parameters_like(components::Array{String, 1}, selected_method, user_input::AbstractString="None")
    filepath = joinpath("database", selected_method, "data_" * selected_method * "_like" * ".csv")
    header = read_line(filepath, 3)
    found_method = search_database_like(components, selected_method)
    found_parameters = Dict{Set{String}, Dict{String, Any}}() 
    for component in keys(found_method)
        #= found_parameters[component] = Dict() =#
        found_parameters[component] = Dict(zip(header, tryparse_fallback(read_line(filepath, found_method[component][selected_method]))))
    end
    return found_parameters
end

function retrieve_parameters_unlike(components::Array{String, 1}, selected_method, user_input::AbstractString="None")
    filepath = joinpath("database", selected_method, "data_" * selected_method * "_unlike" * ".csv")
    header = read_line(filepath, 3)
    found_method = search_database_unlike(components, selected_method)
    found_parameters = Dict{Set{String}, Dict{String, Any}}() 
    for pair in keys(found_method)
        found_parameters[pair] = Dict(zip(header, tryparse_fallback(read_line(filepath, found_method[pair][selected_method]))))
    end
    return found_parameters
end

function retrieve_parameters(components::Array{String, 1}, selected_method, user_input::AbstractString="None")
    parameters_like = retrieve_parameters_like(components, selected_method)
    parameters_unlike = retrieve_parameters_unlike(components, selected_method)
    parameters_assoc = retrieve_parameters_assoc(components, selected_method)

    pairs = vcat([Set([i]) for i in components], [Set(i) for i in collect(Combinatorics.combinations(components, 2))])
    parameters = Dict{Set{String}, Any}()
    for pair in pairs
        parameters[pair] = Dict()
        if haskey(parameters_like, pair)
            merge!(parameters[pair], parameters_like[pair])
        end
        if haskey(parameters_unlike, pair)
            merge!(parameters[pair], parameters_unlike[pair])
        end
        if haskey(parameters_assoc, pair)
            merge!(parameters[pair], parameters_assoc[pair])
        end
    end
    return parameters
end
