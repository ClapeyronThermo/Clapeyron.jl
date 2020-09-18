function lookup(component::String, method="None")
    methods = ["ogSAFT", "PCSAFT", "sPCSAFT", "SAFTVRSW", "softSAFT", "SAFTVRMie"]
    found_methods = Array{String, 1}()
    found_methods_row = Array{Int32, 1}()

    for i in 1:length(methods)
        count = 1
        f = CSV.Rows("database/data_" * methods[i] * ".csv"; header=3)
        for row in f
            if row.Compound == component
                push!(found_methods, methods[i]) # this can be changed to a Dictionary later
                push!(found_methods_row, count)
                break
            end
            count+=1
        end
    end

    !isempty(found_methods) ? println("Found methods for " * component * " in " * join(found_methods, ", ") * ".") : throw(ErrorException("Parameters for " * component * " do not exist for any of the current methods."))

    if method == "none"
        # just use a method for now if none is provided
        method = found_methods[1]
        row = found_methods_row[1]
    else
        row = method in found_methods ? found_methods_row[findfirst(isequal(method), found_methods)] : throw(ErrorException("Parameters for " * component * " do not exist in method " * method * ". Try using another."))
    end

    println("Using method: " * found_methods[1])


    dF = CSV.read("database/data_" * method * ".csv"; header=3, datarow=row+3, limit=1)

    return method, dF
end
