"""
    sle_solubility(model::CompositeModel, p, T, z; solute)

Calculates the solubility of each component within a solution of the other components, at a given temperature and composition.
Returns a matrix containing the composition of the SLE phase boundary for each component. If `solute` is specified, returns only the solubility of the specified component.

Can only function when solid and liquid models are specified within a CompositeModel.
"""
function sle_solubility(model::CompositeModel,p,T,z;solute=nothing)
    mapping = model.mapping
    if isnothing(solute)
        solute = model.solid.components
    end
    p = p*one(eltype(model))
    T = T*one(eltype(model))
    sol = zeros(length(solute),length(model.components))

    for i in 1:length(solute)
        idx_sol_s = zeros(Bool,length(model.solid.components))
        idx_sol_s[model.solid.components .==solute[i]] .= true
        Tm = model.solid.params.Tm.values[idx_sol_s][1]

        idx_sol_l = zeros(Bool,length(model.fluid.components))
        solute_l = mapping[idx_sol_s][1]
        ν_l = [solute_l[1][i][2] for i in 1:length(solute_l[1])]
        solute_l = [solute_l[1][i][1] for i in 1:length(solute_l[1])]

        idx_sol_l[model.fluid.components .== solute_l] .= true
        idx_solv = zeros(Bool,length(model.fluid.components))
        idx_solv[findfirst(idx_sol_l)] = true

        if T>model.solid.params.Tm.values[idx_sol_s][1]
            error("Temperature above melting point of $(solute[i])")
        end

        solid_r,idx_sol_r = index_reduction(model.solid,idx_sol_s)
        μsol = chemical_potential(solid_r,p,T,[1.])

        zref = Float64.(deepcopy(ν_l))
        zref ./= sum(zref)
        fluid_r,idx_liq_r = index_reduction(model.fluid,idx_sol_l)
        γref = activity_coefficient(fluid_r,p,Tm,zref)
        lnKref = log(prod((γref.*zref).^ν_l))
        μsol[1] += lnKref*Rgas()*T
        x0 = x0_sle_solubility(model,p,T,z,idx_solv,idx_sol_l,ν_l,μsol)
        f!(F,x) = obj_sle_solubility(F,model.fluid,p,T,z[.!(idx_solv)],exp10(x[1]),μsol,idx_sol_l,idx_solv,ν_l)
        results = Solvers.nlsolve(f!,x0,LineSearch(Newton()),NEqOptions(),ForwardDiff.Chunk{1}())
        sol[i,idx_solv] .= exp10(Solvers.x_sol(results)[1])
        sol[i,.!(idx_solv)] = z[.!(idx_solv)]
        sol[i,:] ./= sum(sol[i,:])
    end
    if length(solute) == 1
        return sol[1,:]
    else
        return sol
    end
end

function obj_sle_solubility(F,liquid,p,T,zsolv,solu,μsol,idx_sol_l,idx_solv,ν_l)
    z = zeros(typeof(solu),length(liquid))
    
    z[idx_solv] .= solu
    z[.!(idx_solv)] = zsolv
    z ./= sum(z)
    γliq = activity_coefficient(liquid,p,T,z)
    μliq = Rgas()*T*log.(γliq[idx_sol_l].*z[idx_sol_l])
    μliq = sum(μliq.*ν_l)
    F[1] = μliq - μsol[1]
    return F
end

function x0_sle_solubility(model,p,T,z,idx_solv,idx_sol_l,ν_l,μsol)
    z∞ = zeros(length(z))
    z∞[.!(idx_solv)] = z[.!(idx_solv)]
    z∞[idx_solv] .= 1e-10
    z∞ ./= sum(z∞)
    γ∞ = activity_coefficient(model.fluid,p,T,z∞)[idx_sol_l]
    γ∞ = prod(γ∞.^ν_l)
    x0 = [1/sum(ν_l)*log10(exp(μsol[1]/(Rgas()*T)-log(γ∞)))]
    if x0[1] > 0
        x0 = [log10(0.5)]
    end
    return x0
end