include("csv_parser.jl")

function models_return(selected_model::String = "None")
    # Returns an Array of all models in the database directory if input is "None"
    # If input is any other string, it will return an Array of just that string
    if selected_model == "None"
        return filter(x-> isdir(joinpath(dirname(pathof(OpenSAFT)), "../database", x)), readdir(joinpath(dirname(pathof(OpenSAFT)), "../database")))
    else
        if selected_model in filter(x-> isdir(joinpath(dirname(pathof(OpenSAFT)), "../database",x)), readdir(joinpath(dirname(pathof(OpenSAFT)), "../database")))
            return [selected_model]
        else
            error("Database for selected model $selected_model not found.")
        end
    end
end

function createfilepath(selected_model, database_type; variant="None")
    p1 = splitpath(pathof(OpenSAFT))
    p2 = joinpath(p1[1:end-2]...,"database")

    if variant=="None"
        p3 = ("data_$(selected_model)_$(database_type).csv",)
    else
        p3 = (variant,"data_$(selected_model)_$(variant)_$(database_type).csv")
    end
        res =  joinpath(p2, selected_model,p3...)
    #println(res)
    return res
end

function customdatabase_check(selected_model, database_type, customdatabase; variant="None")
    if !isfile(customdatabase)
        error("Custom database $customdatabase not found.")
    elseif selected_model == "None"
        error("Method has to be specified for custom databases.")
    end
    database_filepath = createfilepath(selected_model, database_type; variant=variant)
    database_header = parseline(database_filepath, 3)
    customdatabase_header = parseline(customdatabase, 3)
    if Set(database_header) != Set(customdatabase_header)
        error("Custom database $customdatabase header not consistent with format in $selected_model (variant: $variant) for type $database_type.")
    end
end


function searchdatabase_like(components::Array, selected_model="None"; customdatabase="None", variant="None")
    # Returns a dictionary of components containing another dictionary
    # where the key is a model in which its database contains that parameter,
    # with the found line number as value
    if customdatabase != "None"
        customdatabase_check(selected_model, "like", customdatabase, variant=variant)
    end
    models = models_return(selected_model)
    found_models = Dict(Set([component]) => Dict{String, Int64}() for component in components)
    for component in components
        for model in models
            filepath = customdatabase == "None" ? createfilepath(model, "like"; variant=variant) : customdatabase
            if isfile(filepath)
                found_lines = findmatches(filepath, component, "species"; header_row=3)
                if length(found_lines) > 1
                    # OpenSAFT takes entry with the highest line number if there are duplicates.
                    #= println("The species $component is not unique in database for model $model. Selecting most recent entry.") =#
                elseif isempty(found_lines)
                    if selected_model != "None"
                        error("The species $component cannot be found in database for model $model.")
                    end
                    continue
                end
                found_models[Set([component])][model] = found_lines[end]
            end
        end
    end
    return found_models
end

function searchdatabase_unlike(components::Array, selected_model="None"; customdatabase="None", variant="None")
    # Returns a dictionary of pairs (no same pair) containing another dictionary
    # where the key is a model in which its database contains that parameter,
    # with the found line number as value
    if customdatabase != "None"
        customdatabase_check(selected_model, "unlike", customdatabase, variant=variant)
    end
    models = models_return(selected_model)
    pairs = [Set(i) for i in collect(combinations(components, 2))]
    found_models = Dict{Set{String}, Dict{String, Int64}}()
    for pair in pairs
        for model in models
            filepath = customdatabase == "None" ? createfilepath(model, "unlike"; variant=variant) : customdatabase
            if isfile(filepath)
                components = [i for i in pair]
                found_lines = findmatches_pair(filepath, components[1], components[2], "species1", "species2"; header_row=3)
                if length(found_lines) > 1
                    # OpenSAFT takes entry with the highest line number if there are duplicates.
                    #= println("The pair (" * join([i for i in pair], ", ") * ") is not unique in database for model " * model * ". Selecting most recent entry.") =#
                end
                if !isempty(found_lines)
                    found_models[pair] = Dict()
                    found_models[pair][model] = found_lines[end]
                end
            end
        end
    end
    return found_models
end

function searchdatabase_assoc(components::Array, selected_model="None"; customdatabase="None", variant="None")
    # Returns a dictionary of pairs (same pair allowed) containing another dictionary
    # where the key is a model in which its database contains that parameter,
    # with the found line number as value
    if customdatabase != "None"
        customdatabase_check(selected_model, "assoc", customdatabase, variant=variant)
    end
    models = models_return(selected_model)
    pairs = vcat([Set([i]) for i in components], [Set(i) for i in collect(Combinatorics.combinations(components, 2))])
    found_models = Dict{Set{String}, Dict{String, Array{Int64,1}}}()
    for pair in pairs
        for model in models
            filepath = customdatabase == "None" ? createfilepath(model, "assoc"; variant=variant) : customdatabase
            if isfile(filepath)
                if length(pair) == 1
                    component = [i for i in pair]
                    components = [component component]
                else
                    components = [i for i in pair]
                end
                found_lines = findmatches_pair(filepath, components[1], components[2], "species1", "species2"; header_row=3)
                if length(found_lines) > 0
                    if !haskey(found_models, pair)
                        found_models[pair] = Dict()
                    end
                    found_models[pair][model] = found_lines
                end
            end
        end
    end
    return found_models
end

function searchdatabase_ideal(components::Array, selected_ideal_model="None"; customdatabase="None", variant="None")
    # Returns a dictionary of components containing another dictionary
    # where the key is a model in which its database contains that parameter,
    # with the found line number as value
    if customdatabase != "None"
        customdatabase_check(selected_ideal_model, customdatabase, variant=variant)
    end
    models = models_return(selected_ideal_model)
    found_models = Dict(Set([component]) => Dict{String, Int64}() for component in components)
    for component in components
        for model in models
            filepath = customdatabase == "None" ? createfilepath(model, "ideal"; variant=variant) : customdatabase
            if isfile(filepath)
                found_lines = findmatches(filepath, component, "species"; header_row=3)
                if length(found_lines) > 1
                    # OpenSAFT takes entry with the highest line number if there are duplicates.
                    #= println("The species $component is not unique in database for model $model. Selecting most recent entry.") =#
                elseif isempty(found_lines)
                    if selected_ideal_model != "None"
                        error("The species $component cannot be found in database for model $model.")
                    end
                    continue
                end
                found_models[Set([component])][model] = found_lines[end]
            end
        end
    end
    return found_models
end

function retrieveparams_like(components::Array, selected_model; customdatabase="None", variant="None", redirect="None")
    # Returns a dictionary of all columns in the like database for a selected model
    # for each component
    filepath = customdatabase == "None" ? createfilepath(selected_model, "like"; variant=variant) : customdatabase
    if redirect != "None"
        selected_model = redirect
    end
    header = parseline(filepath, 3)
    found_model = searchdatabase_like(components, selected_model; customdatabase=customdatabase, variant=variant)
    found_params = Dict{Set{String}, Dict{String, Any}}()
    for component in keys(found_model)
        found_params[component] = Dict(zip(header, parseline(filepath, found_model[component][selected_model])))
    end
    return found_params
end

function retrieveparams_unlike(components::Array, selected_model; customdatabase="None", variant="None", redirect="None")
    # Returns a dictionary of all columns in the unlike database for a selected model
    # for each pair
    filepath = customdatabase == "None" ? createfilepath(selected_model, "unlike"; variant=variant) : customdatabase
    if redirect != "None"
        selected_model = redirect
    end
    header = parseline(filepath, 3)
    found_model = searchdatabase_unlike(components, selected_model; customdatabase=customdatabase, variant=variant)
    found_params = Dict{Set{String}, Dict{String, Any}}()
    for pair in keys(found_model)
        found_params[pair] = Dict(zip(header, parseline(filepath, found_model[pair][selected_model])))
    end
    return found_params
end

function retrieveparams_assoc(components::Array, selected_model; customdatabase="None", variant="None", redirect="None")
    # Returns a dictionary of all columns in the assoc database for a selected model
    # for each pair (same pair allowed)
    filepath = customdatabase == "None" ? createfilepath(selected_model, "assoc"; variant=variant) : customdatabase
    if redirect != "None"
        selected_model = redirect
    end
    header = parseline(filepath, 3)
    found_model = searchdatabase_assoc(components, selected_model; customdatabase=customdatabase, variant=variant)
    pairs = keys(found_model)
    found_params = Dict{Set{Tuple{Set{String}, String}}, Dict{String, Any}}()
    if !isempty(found_model)
        for pair in pairs
            for line_number in found_model[pair][selected_model]
                retrieved = Dict(zip(header, parseline(filepath, line_number)))
                assoc_pair = Set([(Set([retrieved["species1"]]), retrieved["site1"]), (Set([retrieved["species2"]]), retrieved["site2"])])
                found_params[assoc_pair] = retrieved # replace if already present
            end
        end
        #= # Check if reverse pair exists; if not, make create reverse pair equal to original pair =#
        #= for assoc_pair in keys(found_params) =#
        #=     if length(assoc_pair) == 2 =#
        #=         pair_tuples = collect(assoc_pair) =#
        #=         reverse_assoc_pair = Set([(pair_tuples[1][1], pair_tuples[2][2]), (pair_tuples[2][1], pair_tuples[1][2])]) =#
        #=         if !haskey(found_params, reverse_assoc_pair) =#
        #=             found_params[reverse_assoc_pair] = =#
        #=                 found_params[assoc_pair] =#
        #=         end =#
        #=     end =#
        #= end =#
    end
    return found_params
end

function retrieveparams_ideal(components::Array, selected_ideal_model; customdatabase="None", variant="None", redirect="None")
    # Returns a dictionary of all columns in the like database for a selected model
    # for each component
    filepath = customdatabase == "None" ? createfilepath(selected_ideal_model, "ideal"; variant=variant) : customdatabase
    if redirect != "None"
        selected_ideal_model = redirect
    end
    header = parseline(filepath, 3)
    found_model = searchdatabase_ideal(components, selected_ideal_model; customdatabase=customdatabase, variant=variant)
    found_params = Dict{Set{String}, Dict{String, Any}}()
    for component in keys(found_model)
        found_params[component] = Dict(zip(header, parseline(filepath, found_model[component][selected_ideal_model])))
    end
    return found_params
end

function retrieveparams(components::Array, selected_model, selected_ideal_model;
        customdatabase_like = "None", variant_like = "None", redirect_like = "None",
        customdatabase_unlike = "None", variant_unlike = "None", redirect_unlike = "None",
        customdatabase_assoc = "None", variant_assoc = "None", redirect_assoc = "None",
        customdatabase_ideal = "None", variant_ideal = "None", redirect_ideal = "None")
    params_like = retrieveparams_like(components, selected_model; customdatabase = customdatabase_like, variant = variant_like, redirect = redirect_like)
    params_unlike = retrieveparams_unlike(components, selected_model; customdatabase = customdatabase_unlike, variant = variant_unlike, redirect = redirect_unlike)
    params_assoc = retrieveparams_assoc(components, selected_model; customdatabase = customdatabase_assoc, variant = variant_assoc, redirect = redirect_assoc)
    if selected_ideal_model != "Basic"
        params_ideal = retrieveparams_ideal(components, selected_ideal_model; customdatabase = customdatabase_ideal, variant = variant_ideal, redirect = redirect_ideal)
    else
        params_ideal = []
    end
    return [params_like, params_unlike, params_assoc, params_ideal]
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
    assoc_pairs = keys(raw_params_assoc)
    like_params_dict = Dict{String, Dict{Set{String}, Float64}}()
    unlike_params_dict = Dict{String, Dict{Set{String}, Float64}}()
    assoc_params_dict = Dict{String, Dict{Set{Tuple{Set{String}, String}}, Float64}}()
    for like_param in like_params
        like_params_dict[like_param] = Dict()
        for component in components
            param_value = raw_params_like[component][like_param]
            if !ismissing(param_value)
                push!(like_params_dict[like_param], component => param_value)
            end
        end
    end
    for unlike_param in unlike_params
        unlike_params_dict[unlike_param] = Dict()
        for pair in pairs
            if haskey(raw_params_unlike[pair], unlike_param)
                param_value = raw_params_unlike[pair][unlike_param]
                if !ismissing(param_value)
                    push!(unlike_params_dict[unlike_param], pair => param_value)
                end
            end
        end
    end
    for assoc_param in assoc_params
        assoc_params_dict[assoc_param] = Dict()
        for assoc_pair in assoc_pairs
            param_value = raw_params_assoc[assoc_pair][assoc_param]
            if !ismissing(param_value)
                push!(assoc_params_dict[assoc_param], assoc_pair => param_value)
            end
        end
    end
    return like_params_dict, unlike_params_dict, assoc_params_dict
end

function filterparams_ideal(raw_params, ideal_params::T) where T<:Array{String,1}
    # Filters the raw parametrs from retrieveparameters into selected headers
    # Returns dictionaries where the keys are the header of the selected columns
    # One dictionary for each like, unlike, and assoc
    raw_params_ideal = raw_params[4]
    components = keys(raw_params_ideal)

    ideal_params_dict = Dict{String, Dict{Set{String}, Float64}}()
    for ideal_param in ideal_params
        ideal_params_dict[ideal_param] = Dict()
        for component in components
            param_value = raw_params_ideal[component][ideal_param]
            if !ismissing(param_value)
                push!(ideal_params_dict[ideal_param], component => param_value)
            end
        end
    end
    return ideal_params_dict
end
