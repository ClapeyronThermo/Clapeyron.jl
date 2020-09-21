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
    found_methods = Dict(component => Dict{String, Int64}() for component in components)
    for component in components
        for method in methods
            filepath = joinpath("database", method, "data_" * method * "_like" * ".csv")
            if isfile(filepath)
                found_lines = find_match(filepath, component, "species"; header_row=3)
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
                    found_methods[component][method] = found_lines[end]
                end
            end
        end
    end
    return found_methods
end

function search_database_unlike(components::Array{String, 1}, selected_method="None", user_input::AbstractString="None")
    methods = methods_return(selected_method)
    pairs = [Set(i) for i in collect(Combinatorics.combinations(components, 2))]
    found_methods = Dict(pair => Dict{String, Int64}() for pair in pairs)
    for pair in pairs
        for method in methods
            filepath = joinpath("database", method, "data_" * method * "_unlike" * ".csv")
            if isfile(filepath)
                components = [i for i in pair]
                found_lines = find_match_pair(filepath, components[1], components[2], "species1", "species2"; header_row=3)
                if length(found_lines) > 1
                    println("The pair (" * join([i for i in pair], ", ") * ") is not unique in database for method " * method * ". Selecting most recent entry.")
                end
                if length(found_lines) != 0
                    found_methods[pair][method] = found_lines[end]
                end
            end
        end
    end
    return found_methods
end

function retrieve_parameters(components::Array{String, 1}, selected_method, user_input::AbstractString="None")
    method = selected_method
    available_methods_like = search_database_like(components, method)
    available_methods_unlike = search_database_unlike(components, method)
    #= avaliable_methods_assoc = search_database_assoc(components) =#

    filepath_like = joinpath("database", method, "data_" * method * "_like" * ".csv")
    filepath_unlike = joinpath("database", method, "data_" * method * "_unlike" * ".csv")
    filepath_assoc = joinpath("database", method, "data_" * method * "_assoc" * ".csv")

    header_like = read_line(filepath_like, 3)
    header_unlike = read_line(filepath_like, 3)
    header_assoc = read_line(filepath_like, 3)

    if haskey(available_methods_like[components[1]], method)
        #= parameters1_like = Dict(zip(header_like, read_line(filepath_like, available_methods_like[components[1]][method]))) =#
        parameters1_like = Dict(zip(header_like, [tryparse(Float64, i) for i in read_line(filepath_like, available_methods_like[components[1]][method])]))
    else
        parameters1_like = Dict()
    end
    if haskey(available_methods_like[components[2]], method)
        #= parameters2_like = Dict(zip(header_like, read_line(filepath_like, available_methods_like[components[2]][method]))) =#
        parameters2_like = Dict(zip(header_like, [tryparse(Float64, i) for i in read_line(filepath_like, available_methods_like[components[2]][method])]))
    else
        parameters1_like = Dict()
    end
    if haskey(available_methods_unlike, method)
        #= parameters_unlike = Dict(zip(header_unlike, read_line(filepath_unlike, available_methods_unlike[Set(components)][method]))) =#
        parameters_unlike = Dict(zip(header_like, [tryparse(Float64, i) for i in read_line(filepath_unlike, available_methods_unlike[Set(components)][method])]))
    else
        parameters_unlike = Dict()
    end

    return (parameters1_like, parameters2_like, parameters_unlike)
end

#### Functions to try ####

# search_database_like(["methanol", "cyclohexane"])  
# search_database_unlike(["methanol", "cyclohexane"])  
# like1, like2, unlike = retrieve_parameters(["methanol", "cyclohexane"], "PCSAFT") 
