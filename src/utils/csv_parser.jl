#### CSV Reader ####
# These are temporary methods to interpret CSV files
# We will probably revert to using the DelimitedFiles module and use actual CSV
# But for now, we are using semicolon delimitetd DSV

function parseline(filepath::AbstractString, selected_line::Int64)
    # Function to read a line from a CSV
    # Returns a dictionary with Float64 values where possible
    open(filepath) do file
        linecount = 1
        for line in eachline(file)
            if linecount == selected_line
                return tryparse_fallback(split(line, ';'))
            end
            linecount += 1
        end
        error("File only contains " * string(linecount-1) * " lines")
    end
end

function findmatches(filepath::AbstractString, match_string::String, selected_col::Int64 = 1)
    # Function to search for a string in the specified column in a CSV
    # Returns an array of line numbers for where string is found
    found_lines = Array{Int64, 1}()
    open(filepath) do file
        linecount = 1
        for line in eachline(file)
            if match_string == split(line, ';'; limit = selected_col+1)[selected_col]
                push!(found_lines, linecount)
            end
            linecount += 1
        end
        return found_lines
    end
end

function findmatches(filepath::AbstractString, match_string::String, selected_col::String; header_row::Int64 = 1)
    selected_col_numbers = findall(isequal(selected_col), parseline(filepath, header_row))
    if length(selected_col_numbers) == 0
        error("Header " * selected_col * " does not exist.")
    elseif length(selected_col_numbers) > 1
        error("Header " * selected_col * " is not unique")
    end
    return findmatches(filepath, match_string, selected_col_numbers[1])
end

function findmatches_pair(filepath::AbstractString, match_string1::String, match_string2::String,selected_col1::Int64 = 1, selected_col2::Int64 = 2; ordered = false)
    # Function to search for pairs of strings in two specified columns
    # Returns an array of line numbers for where string is found
    found_lines = Array{Int64, 1}()
    match_set = Set([match_string1, match_string2])
    open(filepath) do file
        linecount = 1
        for line in eachline(file)
            row_entries = split(line, ';')
            if ordered
                if match_string1 == row_entries[selected_col1] && match_string2 == row_entries[selected_col2]
                    push!(found_lines, linecount)
                end
            else
                if match_set == Set([row_entries[selected_col1], row_entries[selected_col2]])
                    push!(found_lines, linecount)
                end
            end
            linecount += 1
        end
        return found_lines
    end
end

function findmatches_pair(filepath::AbstractString, match_string1::String, match_string2::String, selected_col1::String, selected_col2::String; header_row::Int64 = 1, ordered = false)
    selected_col_numbers1 = findall(isequal(selected_col1), parseline(filepath, header_row))
    selected_col_numbers2 = findall(isequal(selected_col2), parseline(filepath, header_row))
    if length(selected_col_numbers1) == 0
        error("Header " * selected_col1 * " does not exist.")
    elseif length(selected_col_numbers1) > 1
        error("Header " * selected_col1 * " is not unique")
    end
    if length(selected_col_numbers2) == 0
        error("Header " * selected_col2 * " does not exist.")
    elseif length(selected_col_numbers2) > 1
        error("Header " * selected_col2 * " is not unique")
    end
    return findmatches_pair(filepath, match_string1, match_string2, selected_col_numbers1[1], selected_col_numbers2[1]; ordered = ordered)
end

function tryparse_fallback(list)
    # Wrapper for tryparse to retain string if not parsable to Float, and return missing if empty string
    parsed = [tryparse(Float64, i) for i in list]
    return [list[i] == "" ? missing : parsed[i] == nothing ? list[i] : parsed[i] for i in 1:length(list)]
end

