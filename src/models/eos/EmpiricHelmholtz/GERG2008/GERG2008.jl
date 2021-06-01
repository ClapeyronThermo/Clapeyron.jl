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

#uses extended corresponding states and the propane ancilliary eqs
#naively uses Tc as scaling factor
#ideally we would use a eos-agnostic methodology, the noro-frenkel law could give us some insight
#also,, if an eos agnostic shape factor is implemented, we could use GERG2008 as a provider of shape factors

function x0_sat_pure(model::GERG2008,T)
    Ts = T_scale(model)
    vs = _v_scale(model)*0.001 #remember, vc constants in L/mol
    h = vs*5000.0
    T0 = 369.89*T/Ts
    vl = (1.0/_propaneref_rholsat(T0))*h
    vv = (1.0/_propaneref_rhovsat(T0))*h
    return [log10(vl),log10(vv)]
end

function vcompress_v0(model::GERG2008,x)
    return 1.01*lb_volume(model,x)
end

export GERG2008