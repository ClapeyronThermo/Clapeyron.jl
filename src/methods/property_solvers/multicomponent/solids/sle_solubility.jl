"""
    sle_solubility(model::CompositeModel, p, T, z; solute)

Calculates the solubility of each component within a solution of the other components, at a given temperature, pressure and composition.
Returns a matrix containing the composition of the SLE phase boundary for each component. If `solute` is specified, returns only the solubility of the specified component.

Can only function when solid and fluid models are specified within a CompositeModel.
"""
function sle_solubility(model::CompositeModel,p,T,z;solute=nothing,x0=nothing)
    mapping = model.mapping
    if isnothing(solute)
        solute = model.solid.components
    end
    p = p*one(eltype(model))
    T = T*one(eltype(model))
    sol = zeros(length(solute),length(model.components))
    idxs = convert(Vector{Int},indexin(solute,model.solid.components))
    idx_sol = zeros(Bool,length(model.solid.components))
    idx_sol[idxs] .= true
    for i in 1:length(solute)
        idx_sol_s = zeros(Bool,length(model.solid.components))
        idx_sol_s[model.solid.components .==solute[i]] .= true

        #TODO: express this in terms of melting_temperature
        Tm = model.solid.params.Tm.values[idx_sol_s][1]

        idx_sol_l = zeros(Bool,length(model.fluid.components))
        solute_l = mapping[idx_sol_s][1]
        ν_l = [solute_l[1][i][2] for i in 1:length(solute_l[1])]
        solute_l = [solute_l[1][i][1] for i in 1:length(solute_l[1])]


        for i in solute_l
            idx_sol_l[model.fluid.components .== i] .= true
        end
        idx_solv = zeros(Bool,length(model.fluid.components))
        if length(solute_l) == length(model.fluid)
            idx_solv[findfirst(idx_sol_l)] = true
        else
            idx_solv[.!(idx_sol_l)] .= true
        end

        if T > model.solid.params.Tm.values[idx_sol_s][1]
            error("Temperature above melting point of $(solute[i])")
        end

        solid_r,idx_sol_r = index_reduction(model.solid,idx_sol_s)
        μsol = chemical_potential(solid_r,p,T,[1.])

        zref = 1.0 .* ν_l
        zref ./= sum(zref)
        fluid_r,idx_liq_r = index_reduction(model.fluid,idx_sol_l)
        γref = activity_coefficient(fluid_r,p,Tm,zref)
        Kref = one(eltype(γref))
        for i in 1:length(γref)
            Kref *= (zref[i]*γref[i])^ν_l[i]
        end
        lnKref = log(Kref)
        μsol[1] += lnKref*Rgas()*T
        # println(idx_solv)
        # println(idx_sol_l)



        if isnothing(x0)
            x0 = x0_sle_solubility(model,p,T,z,idx_solv,idx_sol_l,ν_l,μsol)
        end
        μ_ref = reference_chemical_potential(model.fluid,p,T,reference_chemical_potential_type(model.fluid))
        data = (μ_ref,idx_sol_l,idx_solv,μsol[1])
        # println(x0)
        f!(F,x) = obj_sle_solubility(F,model,p,T,z[idx_solv],exp10(x[1]),data,ν_l)
        results = Solvers.nlsolve(f!,x0,LineSearch(Newton()),NEqOptions(),ForwardDiff.Chunk{1}())
        sol[i,.!(idx_solv)] .= exp10(Solvers.x_sol(results)[1])
        sol[i,idx_solv] = z[idx_solv]
        sol[i,:] ./= sum(sol[i,:])
    end
    if length(solute) == 1
        return sol[1,:]
    else
        return sol
    end
end

function obj_sle_solubility(F,model,p,T,zsolv,solu,data,ν_l)
    μ_ref,idx_sol_l,idx_solv,μsol = data
    z = zeros(typeof(solu),length(model.fluid))
    z[.!(idx_solv)] .= solu
    z[idx_solv] .= zsolv
    z ./= sum(z)
    R = Rgas(model.fluid)
    ∑z = sum(z)
    γliq = activity_coefficient(model.fluid,p,T,z/∑z,μ_ref = μ_ref)
    γl = @view(γliq[idx_sol_l])
    zl = @view(z[idx_sol_l])
    μliq = zero(eltype(γliq))
    for i in 1:length(γl)
        xli = zl[i]/∑z
        μliq_i = R*T*log(γl[i]*xli)
        μliq += μliq_i*ν_l[i]
    end
    F[1] = μliq - μsol
    return F
end
function x0_sle_solubility(model,p,T,z,idx_solv,idx_sol_l,ν_l,μsol)
    z∞ = zeros(length(z))
    z∞[idx_solv] .= z[idx_solv]
    z∞[.!(idx_solv)] .= 1e-10
    z∞ ./= sum(z∞)
    γ∞ = activity_coefficient(model.fluid,p,T,z∞)[idx_sol_l]
    γ∞ = prod(γ∞.^ν_l)
    x0 = [1/sum(ν_l)*log10(exp(μsol[1]/(Rgas()*T)-log(γ∞)))]
    # println(x0)

    # if x0[1] > 0
    #     x0 = [log10(0.5)]
    # end
    return x0
end

"""
    sle_solubility_T(model::CompositeModel, z, p=1e5; solute=nothing, x0=nothing)

Calculates the solid-liquid coexistence temperature at composition z. Returns the temperature of the most stable solid complex.

Can only function when solid and fluid models are specified within a CompositeModel.
"""
function sle_solubility_T(model::CompositeModel,z,p=1e5;solute=nothing,x0=nothing)
    mapping = model.mapping
    if isnothing(solute)
        solute = model.solid.components
    end
    p = p*one(eltype(model))
    Tsol = zeros(length(solute),length(model.components))
    idxs = convert(Vector{Int},indexin(solute,model.solid.components))
    idx_sol = zeros(Bool,length(model.solid.components))
    idx_sol[idxs] .= true
    for i in 1:length(solute)
        idx_sol_s = zeros(Bool,length(model.solid.components))
        idx_sol_s[model.solid.components .==solute[i]] .= true

        #TODO: express this in terms of melting_temperature
        Tm = model.solid.params.Tm.values[idx_sol_s][1]

        idx_sol_l = zeros(Bool,length(model.fluid.components))
        solute_l = mapping[idx_sol_s][1]
        ν_l = [solute_l[1][i][2] for i in 1:length(solute_l[1])]
        solute_l = [solute_l[1][i][1] for i in 1:length(solute_l[1])]

        for i in solute_l
            idx_sol_l[model.fluid.components .== i] .= true
        end
        idx_solv = zeros(Bool,length(model.fluid.components))
        if length(solute_l) == length(model.fluid)
            idx_solv[findlast(idx_sol_l)] = true
        else
            idx_solv[.!(idx_sol_l)] .= true
        end

        solid_r,idx_sol_r = index_reduction(model.solid,idx_sol_s)

        zref = 1.0 .* ν_l
        zref ./= sum(zref)
        fluid_r,idx_liq_r = index_reduction(model.fluid,idx_sol_l)

        γref = activity_coefficient(fluid_r,p,Tm,zref)
        Kref = one(eltype(γref))
        for i in 1:length(γref)
            Kref *= (zref[i]*γref[i])^ν_l[i]
        end
        lnKref = log(Kref)


        if isnothing(x0)
            x0 = model.solid.params.Tm.values[idx_sol_s]
        end
        

        data = (idx_sol_l,idx_solv,solid_r,lnKref)
        # println(x0)
        f!(F,x) = obj_sle_solubility_T(F,model,p,x[1],z,data,ν_l)
        results = Solvers.nlsolve(f!,x0,LineSearch(Newton()),NEqOptions(),ForwardDiff.Chunk{1}())
        Tsol[i] = Solvers.x_sol(results)[1]
    end
    return maximum(Tsol)
end

function obj_sle_solubility_T(F,model,p,T,z,data,ν_l)
    idx_sol_l,idx_solv,solid_r,lnKref = data
    μsol = chemical_potential(solid_r,p,T,[1.])[1]
    μsol += lnKref*Rgas(model.fluid)*T

    μ_ref = reference_chemical_potential(model.fluid,p,T,reference_chemical_potential_type(model.fluid))

    R = Rgas(model.fluid)
    ∑z = sum(z)
    γliq = activity_coefficient(model.fluid,p,T,z)
    γl = @view(γliq[idx_sol_l])
    zl = @view(z[idx_sol_l])
    μliq = zero(eltype(γliq))
    for i in 1:length(γl)
        xli = zl[i]/∑z
        μliq_i = R*T*log(γl[i]*xli)
        μliq += μliq_i*ν_l[i]
    end
    F[1] = μliq - μsol
    return F
end