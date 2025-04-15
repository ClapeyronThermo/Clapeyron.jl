abstract type DHModel <: IonModel end

struct DH{ϵ} <: DHModel
    components::Array{String,1}
    RSPmodel::ϵ
    references::Array{String,1}
end

"""
    DH(solvents::Array{String,1},
        ions::Array{String,1};
        RSPmodel = ConstRSP,
        userlocations = String[],
        RSPmodel_userlocations = String[],
        verbose = false)

## Input models
- `RSPmodel`: Relative Static Permittivity Model

## Description
This function is used to create a Debye-Hückel model. The Debye-Hückel term gives the excess Helmholtz energy to account for the electrostatic interactions between ions in solution.

## References
1. Debye, P., Huckel, E. (1923). Phys. Z. 24, 185.
"""
DH

export DH
function DH(solvents,ions; RSPmodel=ConstRSP, userlocations=String[], RSPmodel_userlocations=String[], verbose=false)

    components = deepcopy(ions)
    prepend!(components,solvents)

    references = String[]
    init_RSPmodel = @initmodel RSPmodel(solvents,ions,userlocations = RSPmodel_userlocations, verbose = verbose)

    model = DH(components, init_RSPmodel,references)
    return model
end

function a_res(model::DHModel, V, T, z, iondata)
    return @f(a_dh,iondata)
end

function a_dh(ionmodel::DHModel, V, T, z, iondata)
    Z, σ, ϵ_r = iondata
    return a_dh(V, T, z, Z, σ, ϵ_r)
end

function a_dh(V, T, z, Z, σ, ϵ_r)
    nc = length(Z)
    ∑z = sum(z)
    
    s = e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)
    κ = debye_length(V,T,z,ϵ_r,Z,∑z)

    if iszero(κ)
        return zero(κ)
    end
    #ρ = N_A*sum(z)/V
    #κ = sqrt(4π*s*ρ*sum(z[i]*Z[i]*Z[i] for i ∈ 1:nc)/∑z)
    res = zero(κ)
    for i ∈ @iions
        yi,Zi = σ[i]*κ,Z[i]
        yip1 = yi + 1
        χi = 3/(yi*yi*yi)*(3/2+log1p(yi)-2*yip1+1/2*yip1*yip1)
        res +=z[i]*Zi*Zi*χi
    s = e_c^2/(4π*ϵ_0*ϵ_r*k_B*T)
    I = sum(z[i]*Z[i]*Z[i] for i ∈ model.icomponents)
    κ = Solvers.strong_zero(I) do ii
        sqrt(4π*s*N_A/V)*sqrt(ii)
    end
    #κ = sqrt(4π*s*N_A/V)*sqrt(sum(z[i]*Z[i]*Z[i] for i ∈ model.icomponents))
    #iszero(primalval(κ)) && return zero(κ)
    res = zero(Base.promote_eltype(model,V,T,z))
    for i in model.icomponents
        Zi = Z[i]
        if Z[i] != 0 && !iszero(primalval(z[i]))
            yi = σ[i]*κ
            yip1 = yi + 1
            χi = 3/(yi*yi*yi)*(1.5 + log1p(yi) - 2*yip1 + 0.5*yip1*yip1)
            res +=z[i]*Zi*Zi*χi
        end
    end
    return -1/3*s*κ*res/∑z
    #y = σ*κ
    #χ = @. 3/y^3*(3/2+log1p(y)-2*(1+y)+1/2*(1+y)^2)
    # return -1/3*s*κ*sum(z[i]*Z[i]^2*χ[i] for i ∈ iions)/∑z
end