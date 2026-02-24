struct SaltParam <: ClapeyronParam
    explicit_solvent::Bool
    explicit_components::Vector{String} #neutrals + ions
    implicit_components::Vector{String} #neutrals + salts (ion)
    isalts::Vector{Int} #Indices of the salts
    mat::Matrix{Float64} #used to calculate z(ion) -> z(salt)
    F::LU{Float64, Matrix{Float64}, Vector{Int64}} #used to calculate z(salt) -> z(ion)
end

function explicit_salt_param(comps,salts,Z)
    explicit_solvent = true
    explicit_components = comps
    nions = length(Z)
    implicit_components = Vector{String}(undef,nions - 1)
    mat = zeros(nions,nions)
    isalts = Int[]
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
        push!(isalts,i)
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
    if nneutral < nions
        rr[end] .= Z
    end
    return SaltParam(explicit_solvent,explicit_components,implicit_components,isalts,mat,lu(mat))
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
    F = m.mat
    res = F*z
    resize!(res,length(res) - 1)
    return res
end

function to_ion(m,z::AbstractVector)
    if length(m.isalts) == 0
        zz = similar(z)
        zz .= z
    else
        zz = vcat(z,zero(eltype(z)))
    end
    ldiv!(m.F,zz)
    return zz
end

function to_salt(m,result::FlashResult)
    n_ions = sum(result.fractions) #n mols of ions
    z_bulk = sum(b[i]*x[i] for (b,x) in zip(result.fractions,result.compositions))
    n_salts = sum(salt_compositions(m,z_bulk))

    new_comps = map(Base.Fix1(salt_compositions,m),result.compositions)
    for i in new_comps
        xi = new_comps[i]
        xi ./= sum(xi)
    end
    new_fracs = result.fractions .* n_salts ./ n_ions
    return FlashResult(new_comps,new_fracs,result.volumes,result.data)
end

function to_ion(m,result::FlashResult)
    n_salts = sum(result.fractions) #n mols of ions
    z_bulk = sum(b[i]*x[i] for (b,x) in zip(result.fractions,result.compositions))
    n_ions = sum(ion_compositions(m,z_bulk))

    new_comps = map(Base.Fix1(ion_compositions,m),result.compositions)
    for i in new_comps
        xi = new_comps[i]
        xi ./= sum(xi)
    end
    new_fracs = result.fractions .* n_ions ./ n_salts
    return FlashResult(new_comps,new_fracs,result.volumes,result.data)
end


salt_compositions(m::SaltParam,x) = to_salt(m,x)
ion_compositions(m::SaltParam,w) = to_ion(m,w)


#from a vector of salts, split a SaltParam, returns a salt param and the ion indices
function IS_each_split_model(salt::SaltParam,I_salt)
    nions = length(salt.explicit_components)
    m = salt.mat
    I_ion_bool = Vector{Bool}(undef,nions)
    rr = eachrow(m)
    for i in 1:length(rr)-1
        if in(i,I_salt)
            rri = rr[i]
            for k in 1:nions
                I_ion_bool[k] = !iszero(rri[k])
            end
        end
    end
    I_ion_int = findall(I_ion_bool)
    if length(I_ion_int) == length(I_salt)
        mm = m[I_salt,I_ion_int]
    else
        I_salt_plus_charge = vcat(I_salt,nions)
        mm = m[I_salt_plus_charge,I_ion_int]
    end

    isalts = Int[]
    for i in 1:length(salt.isalts)
        si = salt.isalts[i]
        if salt.isalts[i] in I_salt
            push!(isalts,si)
        end
    end

    split_salt = SaltParam(false,salt.explicit_components[I_ion_int],salt.implicit_components[I_salt],isalts,mm,lu(mm))
    return split_salt,I_ion_int
end

export SaltParam