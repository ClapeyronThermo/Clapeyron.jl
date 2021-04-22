include("consts.jl")
include("model.jl")
include("equations.jl")

idealmodel(model::GERG2008) = model
mw(model::GERG2008) = model.Mw

function Base.show(io::IO, mime::MIME"text/plain",sp::GERG2008)
    ln = length(sp.components)
    println(io,ln,"-element GERG008 model, with compounds:")
    for i in 1:ln-1
        println(io," n",i," : ",sp.components[i])
    end
    print(io," n",ln," : ",sp.components[ln])
end

function Base.show(io::IO,sp::GERG2008)
    print(io,"GERG2008(")
    i = 0
    for name in model.components
        i += 1
        if i > 1
            println(io)
        end
        if miss == false
            print(io," ",name)
        else
            print(io," ",name)
        end
    end
    print(io,")")
end

export GERG2008