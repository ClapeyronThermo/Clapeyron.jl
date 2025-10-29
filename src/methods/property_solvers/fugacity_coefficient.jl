# Function to compute fugacity coefficient
function lnϕ(model::EoSModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))


    RT = Rgas(model)*T
    logZ = log(p*vol/RT/sum(z))

    if cache !== nothing
        result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
        nc = length(lnϕ)
        if nc == 1
            lnϕ[1] = VT_lnϕ_pure(model,vol/sum(z),T,p)
        else
            aux .= 0
            aux[1] = vol
            aux[2:nc+1] = z
            gconfig = Solvers._GradientConfig(hconfig)
            F_res(model, V, T, z) = eos_res(model, V, T, z)
            fun(aux) = F_res(model, aux[1], T, @view(aux[2:nc+1]))
            _result = ForwardDiff.gradient!(result, fun, aux, gconfig, Val{false}())
            dresult = DiffResults.gradient(_result)
            μ_res = @view dresult[2:nc+1]
            lnϕ .= μ_res ./ RT .- logZ
        end
    else
        μ_res = VT_chemical_potential_res(model, vol, T, z)
        lnϕ = μ_res/RT .- logZ
    end
    return lnϕ, vol
end

function lnϕ(model::IdealModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    lnϕ = FillArrays.Zeros(length(z))
    return lnϕ, vol
end

function lnϕ!(cache::Tuple, model::EoSModel, p, T, z=SA[1.];
            phase=:unknown,
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    if model isa IdealModel
        lnϕ .= 0
        return lnϕ,vol
    end
    return lnϕ(model,p,T,z,cache;vol)
end

function lnϕ!(lnϕ::AbstractVector, model::EoSModel, p, T, z=SA[1.],cache = nothing;
            phase=:unknown, 
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    if model isa IdealModel
        lnϕ .= 0
        return lnϕ,vol
    end

    if cache != nothing
        lnϕ(model,p,T,z,cache;vol)
        result,aux,lnϕw,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
        vol = aux[1]
        lnϕ .= lnϕw
    else
        RT = Rgas(model)*T
        μ_res = VT_chemical_potential_res!(lnϕ,model, vol, T, z)
        Z = p*vol/RT/sum(z)
        lnϕ .= μ_res ./ RT .- log(Z)
    end
    return lnϕ, vol
end

function VT_lnϕ_pure(model,V,T,p = pressure(model,V,T))
    RT = Rgas(model)*T
    p_res = p - RT/V
    μ_res = eos_res(model,V,T) + p_res*V
    Z = p*V/RT
    return μ_res/RT - log(Z)
end

function ∂lnϕ_cache(model::EoSModel, p, T, z, ::Val{B}) where B
    TT = Base.promote_eltype(model,p,T,z)
    lnϕ = zeros(TT,length(model))
    aux = zeros(TT,length(model) + 1 + B)
    ∂lnϕ∂n = lnϕ * transpose(lnϕ)
    result = DiffResults.HessianResult(aux)
    ∂lnϕ∂n = lnϕ * transpose(lnϕ)
    ∂lnϕ∂P = similar(lnϕ)
    ∂P∂n = similar(lnϕ)
    hconfig = ForwardDiff.HessianConfig(nothing,result,aux)
    if B
        ∂lnϕ∂T = similar(lnϕ)
    else
        ∂lnϕ∂T = lnϕ
    end
    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig
end

# Function to compute fugacity coefficient and its pressure and composition derivatives
function ∂lnϕ∂n∂P(model::EoSModel, p, T, z=SA[1.], cache = ∂lnϕ_cache(model,p,T,z,Val{false}());
            phase=:unknown, 
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))

    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
    RT = Rgas(model)*T
    V = vol
    n = sum(z)
    Z = p*V/RT/n

    F_res(model, V, T, z) = eos_res(model, V, T, z) / RT
    fun(aux) = F_res(model, aux[1], T, @view(aux[2:(ncomponents+1)]))

    ncomponents = length(z)
    aux[1] = V
    aux[2:end] = z
    result = ForwardDiff.hessian!(result, fun, aux, hconfig, Val{false}())

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
    lnϕ .= ∂F∂n .- log(Z)
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
function ∂lnϕ∂n∂P∂T(model::EoSModel, p, T, z=SA[1.],cache = ∂lnϕ_cache(model,p,T,z,Val{true}());
            phase=:unknown, 
            vol0=nothing,
            threaded = true,
            vol = volume(model,p,T,z;phase,vol0,threaded))
    
    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
    RT = Rgas(model)*T
    V = vol
    n = sum(z)
    Z = p*V/RT/n

    ncomponents = length(z)
    aux[1] = V
    aux[2] = T
    aux[3:end] .= z
    F_res(model, V, T, z) = eos_res(model, V, T, z) / RT
    fun(aux) = F_res(model, aux[1], aux[2], @view(aux[3:(ncomponents+2)]))
    result = ForwardDiff.hessian!(result, fun, aux, hconfig, Val{false}())

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
