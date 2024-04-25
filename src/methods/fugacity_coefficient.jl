# Function to compute fugacity coefficient
function lnϕ(model::EoSModel, p, T, z=SA[1.]; phase=:unknown, vol0=nothing)
    RT = R̄*T
    # vol0 === nothing && (vol0 = x0_volume(model, p, T, z, phase = phase))
    # vol = _volume_compress(model,p,T,z,vol0)
    vol = volume(model, p, T, z, phase=phase, vol0=vol0)
    μ_res = VT_chemical_potential_res(model, vol, T, z)
    Z = p*vol/RT/sum(z)
    lnϕ = μ_res/RT .- log(Z)
    return lnϕ, vol
end

function lnϕ!(lnϕ, model::EoSModel, p, T, z=SA[1.]; phase=:unknown, vol0=nothing)
    RT = R̄*T
    # vol0 === nothing && (vol0 = x0_volume(model, p, T, z, phase = phase))
    # vol = _volume_compress(model,p,T,z,vol0)
    vol = volume(model, p, T, z, phase=phase, vol0=vol0)
    μ_res = VT_chemical_potential_res!(lnϕ,model, vol, T, z)
    Z = p*vol/RT/sum(z)
    lnϕ .= μ_res ./ RT .- log(Z)
    return lnϕ, vol
end


function ∂lnϕ_cache(model::EoSModel, p, T, z, dt::Val{B}) where B 
    V = p
    if B
        aux = vcat(V,T,z)
    else
        aux = vcat(V,z)
    end
    lnϕ = zeros(@f(Base.promote_eltype),length(model))
    aux = zeros(@f(Base.promote_eltype),length(model) + 1)
    ∂lnϕ∂n = lnϕ * transpose(lnϕ)
    result = DiffResults.HessianResult(aux)
    ∂lnϕ∂n = lnϕ * transpose(lnϕ)
    ∂lnϕ∂P = similar(lnϕ)
    ∂P∂n = similar(lnϕ)
    if B
        ∂lnϕ∂T = similar(lnϕ)  
    else
        ∂lnϕ∂T = lnϕ
    end
    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T
end

# Function to compute fugacity coefficient and its pressure and composition derivatives
function ∂lnϕ∂n∂P(model::EoSModel, p, T, z=SA[1.],cache = ∂lnϕ_cache(model,p,T,z,Val{false}()); phase=:unknown, vol0=nothing)
    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T = cache
    RT = R̄*T
    V = volume(model, p, T, z, phase=phase, vol0=vol0)
    n = sum(z)
    Z = p*V/RT/n

    F_res(model, V, T, z) = eos_res(model, V, T, z) / RT
    fun(aux) = F_res(model, aux[1], T, @view(aux[2:(ncomponents+1)]))

    ncomponents = length(z)
    aux[1] = V
    aux[2:end] = z
    result = DiffResults.HessianResult(aux)
    result = ForwardDiff.hessian!(result, fun, aux)

    F = DiffResults.value(result)
    ∂F = DiffResults.gradient(result)
    ∂2F = DiffResults.hessian(result)

    ∂F∂V = ∂F[1]
    ∂F∂n = @view ∂F[2:(ncomponents+1)]

    ∂2F∂V2 = ∂2F[1, 1]
    ∂2F∂n2 = @view ∂2F[2:(ncomponents+1), 2:(ncomponents+1)]
    ∂2F∂n∂V = @view ∂2F[1, 2:(ncomponents+1)]

    ∂P∂V = -RT*∂2F∂V2 - n*RT/V^2
    ∂P∂n .= - RT .* ∂2F∂n∂V .+ RT ./ V
    lnϕ = ∂F∂n .- log(Z)
    for i in 1:ncomponents
        ∂P∂ni = ∂P∂n[i]
        ∂V∂ni = - ∂P∂ni/∂P∂V
        ∂lnϕ∂P[i] = ∂V∂ni/RT - 1/p
        for j in 1:ncomponents
            ∂lnϕ∂n[i,j] = ∂2F∂n2[i,j] + 1/n + (∂P∂ni * ∂P∂n[j])/∂P∂V/RT
        end
    end
    return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, V
end


# Function to compute fugacity coefficient and its temperature, pressure and composition derivatives
function ∂lnϕ∂n∂P∂T(model::EoSModel, p, T, z=SA[1.],cache = ∂lnϕ_cache(model,p,T,z,Val{true}()); phase=:unknown, vol0=nothing)
    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T = cache
    RT = R̄*T
    V = volume(model, p, T, z, phase=phase, vol0=vol0)
    n = sum(z)
    Z = p*V/RT/n

    ncomponents = length(z)
    aux[1] = V
    aux[2] = T
    aux[2:end] = z
    F_res(model, V, T, z) = eos_res(model, V, T, z) / R̄ / T
    fun(aux) = F_res(model, aux[1], aux[2], @view(aux[3:(ncomponents+2)]))

    result = DiffResults.HessianResult(aux)
    result = ForwardDiff.hessian!(result, fun, aux)

    F = DiffResults.value(result)
    ∂F = DiffResults.gradient(result)
    ∂2F = DiffResults.hessian(result)

    ∂F∂V = ∂F[1]
    ∂F∂T = ∂F[2]
    ∂F∂n = @view ∂F[3:(ncomponents+2)]

    ∂2F∂V2 = ∂2F[1, 1]
    ∂2F∂n2 = @view ∂2F[3:(ncomponents+2), 3:(ncomponents+2)]
    ∂2F∂T2 = ∂2F[2, 2]
    ∂2F∂n∂V = @view ∂2F[1, 3:(ncomponents+2)]
    ∂2F∂n∂T = @view ∂2F[2, 3:(ncomponents+2)]
    ∂2F∂V∂T = ∂2F[1, 2]

    lnϕ .= ∂F∂n .- log(Z)

    ∂P∂V = -RT*∂2F∂V2 - n*RT/V^2
    ∂P∂n .= -RT .* ∂2F∂n∂V .+ RT ./ V
    ∂P∂T = -RT*∂2F∂V∂T + p/T

    for i in 1:ncomponents
        ∂P∂ni = ∂P∂n[i]
        ∂V∂ni = - ∂P∂ni/∂P∂V
        ∂lnϕ∂P[i] = ∂V∂ni/RT - 1/p
        ∂lnϕ∂T[i] = ∂2F∂n∂T[i] + 1/T - (∂V∂ni*∂P∂T)/RT
        for j in 1:ncomponents
            ∂lnϕ∂n[i,j] = ∂2F∂n2[i,j] + 1/n + (∂P∂ni * ∂P∂n[j])/∂P∂V/RT
        end
    end
    return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, ∂lnϕ∂T, V
end
