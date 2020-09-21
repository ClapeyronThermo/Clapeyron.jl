#### CSV Reader ####
# These are temporary methods to interpret CSV files
# We will probably revert to using the DelimitedFiles module and use actual CSV
# But for now, we are using semicolon delimitetd DSV

function read_line(filename::AbstractString, selected_line::Int64)
    # Function to read a line from a CSV
    # Returns an array of strings
    # Parse to floats later with parse(Float64, string)
    if !isfile(filename)
        error("File does not exist")
    end
    open(filename) do file
        linecount = 1
        for line in eachline(file)
            if linecount == selected_line
                return split(line, ';')
            end
            linecount += 1
        end
        error("File only contains " * string(linecount-1) * " lines")
    end
end

function find_match(filename::AbstractString, match_string::String, selected_col::Int64 = 1)
    # Function to search for a string in the specified column in a CSV
    if !isfile(filename)
        error("File does not exist")
    end
    found_lines = Array{Int64, 1}()
    open(filename) do file
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

function find_match(filename::AbstractString, match_string::String, selected_col::String; header_row::Int64 = 1)
    if !isfile(filename)
        error("File does not exist")
    end
    selected_col_numbers = findall(isequal(selected_col), read_line(filename, header_row))
    if length(selected_col_numbers) == 0
        error("Header " * selected_col * " does not exist.")
    elseif length(selected_col_numbers) > 1
        error("Header " * selected_col * " is not unique")
    end
    return find_match(filename, match_string, selected_col_numbers[1])
end

function find_match_pair(filename::AbstractString, match_string1::String, match_string2::String,selected_col1::Int64 = 1, selected_col2::Int64 = 2)
    # Function to search for pairs of strings in two specified columns
    if !isfile(filename)
        error("File does not exist")
    end
    found_lines = Array{Int64, 1}()
    match_set = Set([match_string1, match_string2])
    open(filename) do file
        linecount = 1
        for line in eachline(file)
            row_entries = split(line, ';')
            if match_set == Set([row_entries[selected_col1], row_entries[selected_col2]])
                push!(found_lines, linecount)
            end
            linecount += 1
        end
        return found_lines
    end
end

function find_match_pair(filename::AbstractString, match_string1::String, match_string2::String, selected_col1::String, selected_col2::String; header_row::Int64 = 1)
    if !isfile(filename)
        error("File does not exist")
    end
    selected_col_numbers1 = findall(isequal(selected_col1), read_line(filename, header_row))
    selected_col_numbers2 = findall(isequal(selected_col2), read_line(filename, header_row))
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
    return find_match_pair(filename, match_string1, match_string2, selected_col_numbers1[1], selected_col_numbers2[1])
end
