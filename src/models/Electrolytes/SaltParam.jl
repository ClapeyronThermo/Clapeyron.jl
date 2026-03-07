struct SaltParam <: ClapeyronParam
    explicit_solvent::Bool
    explicit_components::Vector{String} #neutrals + ions
    implicit_components::Vector{String} #neutrals + salts (ion)
    isalts::Vector{Int} #Indices of the salts
    ion_mat::Matrix{Float64} #used to calculate z(ion) -> z(salt)
    salt_mat::Matrix{Float64} #used to calculate z(salt) -> z(ion)
    E::Matrix{Bool} #stoichiometric matrix
end

function explicit_salt_param(comps,salts,Z)
    explicit_solvent = true
    explicit_components = comps
    ncomps = length(Z)
    nneutral = count(iszero,Z)
    nions = ncomps - nneutral
    nsalts = max(nneutral,nneutral + nions - 1)
    implicit_components = Vector{String}(undef,nsalts)
    salt_mat = zeros(ncomps,ncomps)
    E = zeros(Bool,ncomps,nsalts)
    isalts = Int[]
    #we suppose that first there are nneutral neutral components, followed by nions - nneutral ions
    for i in 1:nneutral
        salt_mat[i,i] = 1
        E[i,i] = true
        implicit_components[i] = explicit_components[i]
    end
    salt_matk = eachcol(salt_mat)
    Ek = eachcol(E)
    k = 0
    for i in (nneutral+1):(nsalts)
        k += 1
        salt = salts[k]
        salt_component = first(salt)
        implicit_components[i] = first(salt)
        push!(isalts,i)
        pairings = last(salt)
        salt_mat_i = salt_matk[i]
        ei = Ek[i]
        for ion_vals in pairings
            ion_i,ni = first(ion_vals),last(ion_vals)
            ki = findfirst(isequal(ion_i),comps)
            if !isnothing(ki)
                ei[ki] = true
                salt_mat_i[ki] = ni
            else
                throw(error("cannot find ions in the salt $salt_component"))
            end
        end
    end
    if ncomps != nneutral
        salt_matk[end] .= Z
    end
    ion_mat = inv(lu(salt_mat))
    return SaltParam(explicit_solvent,explicit_components,implicit_components,isalts,ion_mat,salt_mat,E)
end

SaltParam(model::ESElectrolyteModel) = SaltParam(model,nothing)

function SaltParam(model::ESElectrolyteModel,::Nothing)
    explicit_salt_param(component_list(model),auto_binary_salts(model),model.charge)
end

function SaltParam(model::ESElectrolyteModel,salts)
    explicit_salt_param(component_list(model),salts,model.charge)
end

function component_list(m::SaltParam)
    if m.explicit_solvent
        return m.explicit_components
    else
        return m.implicit_components
    end
end

function to_salt(m,z::AbstractVector)
    zz = m.ion_mat*z
    return deleteat!(zz,length(zz))
end

function to_ion(m,z::AbstractVector)
    if length(m.implicit_components) == length(m.explicit_components)
        zz = similar(z)
        zz .= z
        return zz
    else
        S = @view m.salt_mat[:,1:end-1]
        return S*z
    end
end

function to_salt(m,result::FlashResult)
    new_comps = similar(result.compositions)
    new_volumes = similar(result.volumes)
    new_fracs = similar(result.fractions)
    for i in 1:length(new_comps)
        xi = result.compositions[i]
        wi = salt_compositions(m,xi)
        nwi = sum(wi)
        wi ./= nwi
        new_comps[i] = wi
        new_fracs[i] = result.fractions[i]/nwi
        new_volumes[i] = result.volumes[i]/nwi
    end
    return FlashResult(new_comps,new_fracs,new_volumes,result.data)
end

function to_ion(m,result::FlashResult)
    new_comps = similar(result.compositions)
    new_volumes = similar(result.volumes)
    new_fracs = similar(result.fractions)
    for i in 1:length(new_comps)
        xi = result.compositions[i]
        wi = ion_compositions(m,xi)
        nwi = sum(wi)
        wi ./= nwi
        new_comps[i] = wi
        new_fracs[i] = result.fractions[i]/nwi
        new_volumes[i] = result.volumes[i]/nwi
    end
    return FlashResult(new_comps,new_fracs,new_volumes,result.data)
end


salt_compositions(m::SaltParam,x) = to_salt(m,x)
ion_compositions(m::SaltParam,w) = to_ion(m,w)


#from a vector of salts, split a SaltParam, returns a salt param and the ion indices
function IS_each_split_model(salt::SaltParam,I_salt)
    nions = length(salt.explicit_components)
    nsalts = length(salt.implicit_components)
    I_ion_bool = Vector{Bool}(undef,nions)
    I_ion_bool .= false
    E = salt.E
    EE = eachcol(E)
    for i in 1:nions
        if in(i,I_salt)
            Ei = EE[i]
            for k in 1:nions
                I_ion_bool[k] = I_ion_bool[k] | Ei[k]
            end
        end
    end
    
    I_ion_int = findall(I_ion_bool)
    if length(I_ion_int) == length(I_salt)
        salt_mat = salt.salt_mat[I_ion_int,I_salt]
    else
        I_salt_plus_charge = vcat(I_salt,nions)
        salt_mat = salt.salt_mat[I_ion_int,I_salt_plus_charge]
    end
    EE = salt.E[I_ion_int,I_salt]
    isalts = Int[]
    for i in 1:length(salt.isalts)
        si = salt.isalts[i]
        if salt.isalts[i] in I_salt
            push!(isalts,si)
        end
    end
    ion_mat = inv(lu(salt_mat))
    split_salt = SaltParam(false,salt.explicit_components[I_ion_int],salt.implicit_components[I_salt],isalts,ion_mat,salt_mat,EE)
    return split_salt,I_ion_int
end

export SaltParam
