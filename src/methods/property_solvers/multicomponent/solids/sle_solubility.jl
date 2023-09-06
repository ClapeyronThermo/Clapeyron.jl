"""
    sle_solubility(model::CompositeModel, p, T, z; solute)

Calculates the solubility of each component within a solution of the other components, at a given temperature and composition.
Returns a matrix containing the composition of the SLE phase boundary for each component. If `solute` is specified, returns only the solubility of the specified component.

Can only function when solid and liquid models are specified within a CompositeModel.
"""
function sle_solubility(model::CompositeModel,p,T,z;solute=nothing)
    if isnothing(solute)
        solute = model.components
    end
    
    sol = zeros(length(solute),length(model.components))

    for i in 1:length(solute)
        
        idx_sol = zeros(Bool,length(model.components))
        idx_sol[model.components .==solute[i]] .= true

        if T>model.solid.params.Tm.values[idx_sol][1]
            error("Temperature above melting point of $(solute[i])")
        end

        solid_r,idx_sol_r = index_reduction(model.solid,idx_sol)
        μsol = chemical_potential(solid_r,p,T,[1.])
        
        x0 = x0_sle_solubility(model,p,T,z,idx_sol,μsol)
        # println(x0)
        f!(F,x) = obj_sle_solubility(F,model.liquid,p,T,z[.!(idx_sol)],exp10(x[1]),μsol,idx_sol)
        results = Solvers.nlsolve(f!,x0,LineSearch(Newton()),NEqOptions(),ForwardDiff.Chunk{1}())
        sol[i,idx_sol] .= exp10(Solvers.x_sol(results)[1])
        sol[i,.!(idx_sol)] = z[.!(idx_sol)]
        sol[i,:] ./= sum(sol[i,:])
    end
    if length(solute) == 1
        return sol[1,:]
    else
        return sol
    end
end

function obj_sle_solubility(F,liquid,p,T,zsolv,solu,μsol,idx_sol)
    z = zeros(typeof(solu),length(zsolv)+1)
    z[.!(idx_sol)] = zsolv
    z[idx_sol] .= solu
    z ./= sum(z)
    γliq = activity_coefficient(liquid,p,T,z)
    μliq = Rgas()*T*log(γliq[idx_sol][1]*z[idx_sol][1])
    F[1] = μliq - μsol[1]
    return F
end

function x0_sle_solubility(model,p,T,z,idx_sol,μsol)
    z∞ = zeros(length(z))
    z∞[.!(idx_sol)] = z[.!(idx_sol)]
    z∞[idx_sol] .= 1e-10
    z∞ ./= sum(z∞)
    γ∞ = activity_coefficient(model.liquid,p,T,z∞)[idx_sol][1]
    x0 = [log10(exp(μsol[1]/(Rgas()*T)-log(γ∞)))]
    return x0
end