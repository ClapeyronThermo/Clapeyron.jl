function Base.show(io::IO, ::MIME"text/plain", model::NonGCSAFTModel)
    if get(io, :compact, false)
        println(model.modelname, " model with using ideal model ", model.idealmodel.modelname, " with ", model.lengthcomponents, " components:")
        for i in model.icomponents
            print(" ", model.components[i])
            println()
        end
        print("Contains parameters: ")
        firstloop = true
        for fieldname in fieldnames(typeof(model.params))
            firstloop == false && print(", ")
            print(getfield(model.params, fieldname).name)
            firstloop = false
        end
    else
        print(model.modelname, "{", model.idealmodel.modelname, "}(")
        firstloop = true
        for i in model.icomponents
            firstloop == false && print(", ")
            print(model.components[i])
            firstloop = false
        end
        print(")")
    end
end

function Base.show(io::IO, ::MIME"text/plain", model::GCSAFTModel)
    if get(io, :compact, false)
        println(model.modelname, " model with using ideal model ", model.idealmodel.modelname, " with ", model.lengthcomponents, " components:")
        for i in model.icomponents
            print(" ", model.components[i], ": ")
            firstloop = true
            for k in 1:length(model.allcomponentgroups[i])
                firstloop == false && print(", ")
                print(model.allcomponentgroups[i][k], " => ", model.allcomponentngroups[i][k])
                firstloop = false
            end
            println()
        end
        print("Contains parameters: ")
        firstloop = true
        for fieldname in fieldnames(typeof(model.params))
            firstloop == false && print(", ")
            print(getfield(model.params, fieldname).name)
            firstloop = false
        end
    else
        print(model.modelname, "{", model.idealmodel.modelname, "}(")
        firstloop = true
        for i in model.icomponents
            firstloop == false && print(", ")
            print(model.components[i])
            firstloop = false
        end
        print(")")
    end
end
