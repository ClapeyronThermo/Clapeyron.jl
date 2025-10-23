"""
    eutectic_point(model::CompositeModel, p)

Calculates the eutectic point of a binary mixture (at a given pressure).
Returns a tuple containing the eutectic temperature and the composition at the eutectic point.

Can only function when solid and liquid models are specified within a `CompositeModel`.
"""
function eutectic_point(model::CompositeModel,p=1e5; x0=nothing)
    p = p*one(eltype(model))
    binary_component_check(eutectic_point,model)
    solid = solid_model(model)
    fluid = fluid_model(model)
    f!(F,x) = obj_eutectic_point(F,solid,fluid,p,x[1]*200.,FractionVector(x[2]))
    x0 === nothing && (x0 = x0_eutectic_point(model,p))
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
    γliq = activity_coefficient(liquid,p,T,x; phase=:l)
    RT = Rgas()*T
    μliq = @. Rgas()*T*log(γliq*x)
    dμ1 = exp(μsol[1]/RT)/(γliq[1]*x[1])
    dμ2 = exp(μsol[2]/RT)/(γliq[2]*x[2])
    F[1] = dμ1
    F[2] = dμ2
    #F[1:2] = μliq .- μsol
    return F
end

function x0_eutectic_point(model::EoSModel,p)
    pure = split_pure_model(model)
    fus = melting_temperature.(pure,p)
    Tm = first.(fus)
    vs = getindex.(fus,2)
    vl = last.(fus)
    K(modeli,_Vs,_Vl,_T) = -dpdT_saturation(fluid_model(modeli),solid_model(modeli),_Vl,_Vs,_T)*_T*_T/p
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