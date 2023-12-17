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
    f!(F,x) = obj_eutectic_point(F,model.solid,model.fluid,p,x[1]*200.,FractionVector(x[2]))
    x0 = x0_eutectic_point(model)
    # println(x0)
    results = Solvers.nlsolve(f!,x0)
    T = Solvers.x_sol(results)[1]*200
    x = FractionVector(Solvers.x_sol(results)[2])
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

    _Hfus(modeli,_Vs,_Vl,_T,z =SA[1.0]) = VT_enthalpy(fluid_model(model),_Vl,_T,SA[1.0]) - VT_enthalpy(solid_model(model),_Vs,_T,SA[1.0])
    Hfus = _Hfus.(pure,vs,vl,Tm)
    #Hfus correlation
    Tmax = 2/minimum(Tm)
    f(tinv) = 1 - exp(-Hfus[1]/R*(tinv -1/Tm[1])) - exp(-Hfus[2]/R*(tinv -1/Tm[2]))
    prob = Roots.ZeroProblem(f,(zero(Tmax),Tmax))
    T0inv = Roots.solve(prob)
    T0 = 1/T0inv
    x0 = exp(-Hfus[1]/Rgas()*(1/T0-1/Tm[1]))
    return (T0/200.,x0)
end