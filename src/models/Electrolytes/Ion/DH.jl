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
    
    κ = debye_length(V,T,z,ϵ_r,Z)
    res = zero(Base.promote_eltype(κ,σ))
    for i in 1:nc
        Zi = Z[i]
        if Z[i] != 0 && !iszero(primalval(z[i]))
            σi = σ[i]
            yi = σ[i]*κ
            yip1 = yi + 1
            #χi = 3/(yi*yi*yi)*(1.5 + log1p(yi) - 2*yip1 + 0.5*yip1*yip1)
            χi = (log1p(yi) + 0.5*yi*(yi - 2))/(yi*yi*yi)
            res +=z[i]*Zi*Zi*χi
        end
    end
    s = e_c*e_c/(4π*ϵ_0*ϵ_r*k_B*T)
    return -1*s*κ*res/∑z
    #y = σ*κ
    #χ = @. 3/y^3*(3/2+log1p(y)-2*(1+y)+1/2*(1+y)^2)
    # return -1/3*s*κ*sum(z[i]*Z[i]^2*χ[i] for i ∈ iions)/∑z
end