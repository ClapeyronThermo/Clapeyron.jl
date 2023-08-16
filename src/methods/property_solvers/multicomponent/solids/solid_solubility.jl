function solid_solubility(model::CompositeModel,solute::String,p,T,zsolv)
    components = model.components
    idx_sol = solute.==components
    solid_r,idx_sol_r = index_reduction(model.solid,idx_sol)
    μsol = chemical_potential(solid_r,p,T,[1.])
    f!(F,x) = obj_solid_solubility(F,model.liquid,p,T,zsolv,exp10(x[1]),μsol,idx_sol)
    x0 = [-2.]
    results = Solvers.nlsolve(f!,x0,LineSearch(Newton()),NEqOptions(),ForwardDiff.Chunk{1}())
    sol = exp10(Solvers.x_sol(results)[1])
    return sol
end

function obj_solid_solubility(F,liquid,p,T,zsolv,solu,μsol,idx_sol)
    z = vcat(zsolv,solu)
    z ./= sum(z)
    γliq = activity_coefficient(liquid,p,T,z)
    μliq = log(γliq[idx_sol][1]*z[idx_sol][1])
    F[1] = μliq - μsol[1]
    return F
end