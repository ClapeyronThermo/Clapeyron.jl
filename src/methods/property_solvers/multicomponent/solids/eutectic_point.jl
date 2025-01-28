"""
    eutectic_point(model::CompositeModel, p)

Calculates the eutectic point of a binary mixture (at a given pressure).
Returns a tuple containing the eutectic temperature and the composition at the eutectic point.

Can only function when solid and liquid models are specified within a CompositeModel.
"""
function eutectic_point(model::CompositeModel,p=1e5)
    p = p*one(eltype(model))
    if length(model) != 2
        error("Eutectic point only defined for binary systems")
    end
    solid = solid_model(model)
    fluid = fluid_model(model)
    f!(F,x) = obj_eutectic_point(F,solid,fluid,p,x[1]*200.,FractionVector(x[2]))
    x0 = x0_eutectic_point(model,p)
    # println(x0)
    r  = Solvers.nlsolve(f!,x0)
    sol = Solvers.x_sol(r)
    !all(<(r.options.f_abstol),r.info.best_residual) && (sol .= NaN)
    T = sol[1]*200
    x = FractionVector(sol[2])
    return T,x
end

function obj_eutectic_point(F,solid,liquid,p,T,x)
    μsol = chemical_potential(solid,p,T,x)
    γliq = activity_coefficient(liquid,p,T,x)
    μliq = @. Rgas()*T*log(γliq*x)
    F[1:2] = μliq .- μsol
    return F
end

function x0_eutectic_point(model::EoSModel,p)
    pure = split_model(model)
    fus = melting_temperature.(pure,p)
    Tm = first.(fus)
    vs = getindex.(fus,2)
    vl = last.(fus)
    K(modeli,_Vs,_Vl,_T) = -dpdT_pure(fluid_model(modeli),solid_model(modeli),_Vl,_Vs,_T)*_T*_T/p
    Ki = K.(pure,vs,vl,Tm)
    #Clausius-Clapeyron correlation
    Tmax = 2/minimum(Tm)
    f(tinv) = 1 - exp(Ki[1]*(tinv -1/Tm[1])) - exp(Ki[2]*(tinv -1/Tm[2]))
    prob = Roots.ZeroProblem(f,(zero(Tmax),Tmax))
    T0inv = Roots.solve(prob)
    T0 = 1/T0inv
    x0 = exp(Ki[1]*(1/T0-1/Tm[1]))
    return [T0/200.,x0]
end