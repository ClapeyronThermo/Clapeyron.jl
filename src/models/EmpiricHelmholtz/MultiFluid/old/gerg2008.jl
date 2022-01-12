#=
idealmodel(model::GERG2008) = model
molecular_weight(model::GERG2008,z=SA[1.0]) = comp_molecular_weight(mw(model),z)
split_model(model::GERG2008) = simple_split_model(model)
mw(model::GERG2008) = model.Mw
function Base.show(io::IO, mime::MIME"text/plain",sp::GERG2008)
    ln = length(sp.components)
    println(io,"GERG008 model with ",ln," component",ifelse(isone(ln),"","s"),":")
    tick = "\""
    for i in 1:ln-1
        println(io,tick,sp.components[i],tick)
    end
    print(io,tick,sp.components[ln],tick)
end

function Base.show(io::IO,model::GERG2008)
    return eosshow(io,model)
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

function x0_volume_liquid(model::GERG2008,T,z)
    return 1.01*lb_volume(model,z)
end

export GERG2008
=#