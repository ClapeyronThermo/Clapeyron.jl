function LLE(model::ActivityModel,T;v0=nothing)
    if v0 === nothing
        if length(model) == 2
        v0 = [0.25,0.75]
        else
            throw(error("unable to provide an initial point for LLE pressure"))
        end
    end
    
    len = length(v0)

    Fcache = zeros(eltype(v0),len)
    f!(F,z) = Obj_LLE(model, F, T, z[1], z[2])
    r  = Solvers.nlsolve(f!,v0,LineSearch(Newton()))
    sol = Solvers.x_sol(r)
    x = sol[1]
    xx = sol[2]
    return x,xx
end

function Obj_LLE(model::ActivityModel, F, T, x, xx)
    x = Fractions.FractionVector(x)
    xx = Fractions.FractionVector(xx)
    γₐ = activity_coefficient(model,1e-3,T,x)
    γᵦ = activity_coefficient(model,1e-3,T,xx)

    F .= γᵦ.*xx-γₐ.*x
    return F
end

export LLE