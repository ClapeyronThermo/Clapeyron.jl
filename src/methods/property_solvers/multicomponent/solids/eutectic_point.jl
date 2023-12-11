"""
    eutectic_point(model::CompositeModel, p)

Calculates the eutectic point of a binary mixture (at a given pressure).
Returns a tuple containing the eutectic temperature and the composition at the eutectic point.

Can only function when solid and liquid models are specified within a CompositeModel.
"""
function eutectic_point(model::CompositeModel,p=1e5)
    p = p*one(eltype(model))
    if length(model.components) != 2
        error("Eutectic point only defined for binary systems")
    end
    f!(F,x) = obj_eutectic_point(F,model.solid,model.fluid,p,x[1]*200.,FractionVector(x[2]))
    x0 = x0_eutectic_point(model)
    # println(x0)
    results = Solvers.nlsolve(f!,x0)
    return [Solvers.x_sol(results)[1]*200.,FractionVector(Solvers.x_sol(results)[2])]
end

function obj_eutectic_point(F,solid,liquid,p,T,x)
    μsol = chemical_potential(solid,p,T,x)
    γliq = activity_coefficient(liquid,p,T,x)
    μliq = @. Rgas()*T*log(γliq*x)
    F[1:2] = μliq .- μsol
    return F
end

function x0_eutectic_point(model::CompositeModel)
    model_s = model.solid
    Hfus = model_s.params.Hfus.values
    Tm = model_s.params.Tm.values
    f(t) = 1-sum(exp.(-Hfus./Rgas().*(1/t.-1 ./Tm)))
    T0 = Roots.fzero(f,(minimum(Tm)*0.5,minimum(Tm)))
    x0 = exp(-Hfus[1]/Rgas()*(1/T0-1/Tm[1]))
    return [T0/200.,x0]
end