abstract type SolidKsModel <: EoSModel end

struct SolidKsParam <: EoSParam
    Gform::SingleParam{Float64}
    Hform::SingleParam{Float64}
    Tref::SingleParam{Float64}
end

@newmodelsimple SolidKs SolidKsModel SolidKsParam

"""
    SolidKsModel <: EoSModel

    SolidKs(components;
    userlocations = String[],
    verbose::Bool=false)

## Parameters

- `Hfus`: Single Parameter (`Float64`) - Enthalpy of Fusion at 1 bar `[J/mol]`
- `Tm`: Single Parameter (`Float64`) - Melting Temperature `[K]`
- `CpSL`: Single Parameter (`Float64`) - Heat Capacity of the Solid-Liquid Phase Transition `[J/mol/K]`

## Description

Approximation of the excess chemical potential in the solid phase, using enthalpies and gibbs energies of formation:
```
ln(xᵢγᵢ) = -Gformᵢ*T/Trefᵢ - Hformᵢ*(1 - T/Trefᵢ)
```
"""
SolidKs
default_locations(::Type{SolidKs}) = ["solids/formation.csv"]
default_references(::Type{SolidKs}) = String[]
default_ignore_missing_singleparams(::Type{SolidKs}) = ["CpSL"]

function volume_impl(model::SolidKsModel,p,T,z,phase,threaded,vol0)
    _0 = zero(T + first(z))
    return _0/_0
end

sle_T_ref(model::SolidKsModel) = model.params.Tref.values

function chemical_potential_impl(model::SolidKsModel,p,T,z,phase,threaded,vol0)
    Gform = model.params.Gform.values
    Hform = model.params.Hform.values
    Tref = model.params.Tref.values
    return @. -Gform*T/Tref - Hform*(1 - T/Tref)
end

export SolidKs

"""
    sle_solubility(model::CompositeModel, p, T, z; solute)

Calculates the solubility of each component within a solution of the other components, at a given temperature and composition.
Returns a matrix containing the composition of the SLE phase boundary for each component. If `solute` is specified, returns only the solubility of the specified component.

Can only function when solid and fluid models are specified within a CompositeModel.
"""

function sle_solubility(model::CompositeModel{F,S},p,T,z;solute=nothing,x0=nothing) where F <: EoSModel where S <: SolidKsModel
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
        #Tm = model.solid.params.Tm.values[idx_sol_s][1]

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

        #if T > model.solid.params.Tm.values[idx_sol_s][1]
        #    error("Temperature above melting point of $(solute[i])")
        #end

        solid_r,idx_sol_r = index_reduction(model.solid,idx_sol_s)
        μsol = chemical_potential(solid_r,p,T,[1.])

        Tref = sle_T_ref(solid_r)[1]
        zref = 1.0 .* ν_l
        zref ./= sum(zref)
        fluid_r,idx_liq_r = index_reduction(model.fluid,idx_sol_l)
        γref = activity_coefficient(fluid_r,p,Tref,zref)
        Kref = one(eltype(γref))
        for i in 1:length(Kref)
            Kref *= (zref[i]*γref[i])^ν_l[i]
        end
        lnKref = log(Kref)
        #μsol[1] += lnKref*Rgas()*T

        if isnothing(x0)
            x0 = x0_sle_solubility(model,p,T,z,idx_solv,idx_sol_l,ν_l,μsol)
        end
        μ_ref = reference_chemical_potential(model.fluid,p,T,reference_chemical_potential_type(model.fluid))
        data = (μ_ref,idx_sol_l,idx_solv,μsol[1])

        f!(F,x) = obj_sle_solubility(F,model,p,T,z,exp10(x[1]),idx_sol_l,idx_sol_s,idx_solv,ν_l)
        results = Solvers.nlsolve(f!,x0,LineSearch(Newton()),NEqOptions(f_abstol=1e-6,f_reltol=1e-8),ForwardDiff.Chunk{1}())
        sol[i,.!(idx_solv)] .= exp10(Solvers.x_sol(results)[1]).*ν_l
        sol[i,idx_solv] = z[idx_solv]
        sol[i,:] ./= sum(sol[i,:])
    end
    if length(solute) == 1
        return sol[1,:]
    else
        return sol
    end
end

function obj_sle_solubility(F,model::CompositeModel{L,S},p,T,zsolv,solu,idx_sol_l,idx_sol_s,idx_solv,ν_l) where L <: EoSModel where S <: SolidKsModel
    z = zeros(typeof(solu),length(model.fluid))
    z[.!(idx_solv)] .= solu.*ν_l
    z[idx_solv] .= zsolv[idx_solv]
    z ./= sum(z)

    if typeof(model.fluid) <: ESElectrolyteModel
        v = volume(model.fluid,p,T,z)
        μ = VT_chemical_potential_res(model.fluid,v,T,z) .- Rgas()*T*log(v*p/(Rgas()*T*sum(z))) + Rgas()  * T * log.(z)
        
        zref = ones(length(model.fluid))*1e-30
        idx_water = find_water_indx(model.fluid)
        zref[idx_water] = 1.0
        # zref = zeros(length(model.fluid))

        # ineutral = model.fluid.charge .== 0
        # zref[.!(ineutral)] .= 1e-30
        # zref[ineutral] .= zsolv[ineutral]
        zref ./= sum(zref)
        vref = volume(model.fluid,p,T,zref)
        μref = VT_chemical_potential_res(model.fluid,vref,T,zref) .- Rgas()*T*log(vref*p/(Rgas()*T*sum(zref)))

        μliq = (μ - μref)[idx_sol_l]
    else
        pure   = split_pure_model(model.fluid)
        μ_mixt = chemical_potential(model.fluid, p, T, z)
        μ_ref = gibbs_free_energy.(pure, p, T)
        μliq = (μ_mixt - μ_ref)[idx_sol_l]
    end

    solid_r,idx_sol_r = index_reduction(model.solid,idx_sol_s)
    μsol = chemical_potential(solid_r,p,T,[1.])

    μliq = sum(μliq.*ν_l)
    F[1] = (μliq - μsol[1])/(Rgas()*T)
    return F
end

function x0_sle_solubility(model::CompositeModel{L,S},p,T,z,idx_solv,idx_sol_l,ν_l,μsol) where L <: EoSModel where S <: SolidKsModel
    z∞ = zeros(length(z))
    z∞[.!(idx_solv)] = z[.!(idx_solv)]
    z∞[idx_solv] .= 1e-10
    z∞ ./= sum(z∞)
    γ∞ = activity_coefficient(model.fluid,p,T,z∞)[idx_sol_l]
    γ∞ = prod(γ∞.^ν_l)
    x0 = [1/sum(ν_l)*log10(exp(μsol[1]/RT-log(γ∞)))]
    if x0[1] > 0
        x0 = [log10(0.5)]
    end
    return x0
end