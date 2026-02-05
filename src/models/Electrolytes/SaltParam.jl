struct SaltParam <: ClapeyronParam
    explicit_solvent::Bool
    explicit_components::Vector{String} #neutrals + ions
    implicit_components::Vector{String} #neutrals + salts (ion)
    mat::Matrix{Float64} #used to calculate z(ion) -> z(salt)
    F::LU{Float64, Matrix{Float64}, Vector{Int64}} #used to calculate z(salt) -> z(ion)
end

function explicit_salt_param(comps,salts,Z)
    explicit_solvent = true
    explicit_components = comps
    nions = length(Z)
    implicit_components = Vector{String}(undef,nions - 1)
    mat = zeros(nions,nions)
    nneutral = count(iszero,Z)
    #we suppose that first there are nneutral neutral components, followed by nions - nneutral ions
    for i in 1:nneutral
        mat[i,i] = 1
        implicit_components[i] = explicit_components[i]
    end
    rr = eachrow(mat)
    k = 0
    for i in (nneutral+1):(nions-1)
        k += 1
        salt = salts[k]
        salt_component = first(salt)
        implicit_components[i] = first(salt)
        pairings = last(salt)
        ri = rr[i]
        for ion_vals in pairings
            ion_i,ni = first(ion_vals),last(ion_vals)
            ki = findfirst(isequal(ion_i),comps)
            
            if !isnothing(ki)
                ri[ki] = 1/ni
            else
                throw(error("cannot find ions in the salt $salt_component"))
            end
        end
    ∑ri = count(!iszero,ri)
    ri ./=  ∑ri
    end
    rr[end] .= Z
    return SaltParam(explicit_solvent,explicit_components,implicit_components,mat,lu(mat))
end

function SaltParam(model::ESElectrolyteModel)
    explicit_salt_param(component_list(model),auto_binary_salts(model),model.charge)
end

function component_list(m::SaltParam)
    if m.explicit_solvent
        return m.explicit_components
    else
        return m.implicit_components
    end
end

function to_salt(m::SaltParam,z)
    F = m.mat
    res = F*z
    resize!(res,length(res) - 1)
    return res
end

function to_ion(m::SaltParam,z)    
    zz = vcat(z,zero(eltype(z)))
    ldiv!(m.F,zz)
    return zz
end
