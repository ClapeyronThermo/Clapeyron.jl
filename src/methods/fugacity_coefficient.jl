# Function to compute fugacity coefficient
function lnϕ(model::EoSModel, p, T, z=SA[1.],cache = nothing; phase=:unknown, vol0=nothing,threaded = true)
    vol = volume(model, p, T, z, phase=phase, vol0=vol0, threaded=threaded)
    RT = Rgas(model)*T
    logZ = log(p*vol/RT/sum(z))

    if cache != nothing
        result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
        nc = length(lnϕ)
        aux .= 0
        aux[1] = vol
        aux[2:nc+1] = z
        gconf,jconf = hconfig.gradient_config,hconfig.jacobian_config
        seeds = jconf.seeds
        duals = jconf.duals[1]
        gconfig = ForwardDiff.GradientConfig{Nothing,eltype(aux),length(seeds),typeof(duals)}(seeds,duals)
        gresult = ForwardDiff.MutableDiffResult(result.value,(result.derivs[1],))
        F_res(model, V, T, z) = eos_res(model, V, T, z)
        fun(aux) = F_res(model, aux[1], T, @view(aux[2:nc+1]))
        _result = ForwardDiff.gradient!(gresult, fun, aux, gconfig, Val{false}())
        dresult = _result.derivs[1]
        μ_res = @view dresult[2:nc+1]
        lnϕ .= μ_res ./ RT .- logZ
    else
        μ_res = VT_chemical_potential_res(model, vol, T, z)
        Z = p*vol/RT/sum(z)
        lnϕ = μ_res/RT .- log(Z)
    end
    return lnϕ, vol
end

function lnϕ(model::IdealModel, p, T, z=SA[1.],cache = nothing; phase=:unknown, vol0=nothing,threaded = true)
    vol = volume(model, p, T, z, phase=phase, vol0=vol0, threaded=threaded)
    lnϕ = FillArrays.Zeros(length(z))
    return lnϕ, vol
end

function lnϕ!(lnϕ, model::EoSModel, p, T, z=SA[1.]; phase=:unknown, vol0=nothing, threaded = true)
    RT = Rgas(model)*T
    vol = volume(model, p, T, z, phase=phase, vol0=vol0, threaded=threaded)
    μ_res = VT_chemical_potential_res!(lnϕ,model, vol, T, z)
    Z = p*vol/RT/sum(z)
    lnϕ .= μ_res ./ RT .- log(Z)
    return lnϕ, vol
end

function VT_lnϕ_pure(model,V,T,p = pressure(model,V,T))
    RT = Rgas(model)*T
    μ_res = a_res(model,V,T,SA[1.0])
    Z = p*V/RT
    lnϕ = μ_res - log(Z)
end

function ∂lnϕ_cache(model::EoSModel, p, T, z, dt::Val{B}) where B
    V = p
    lnϕ = zeros(@f(Base.promote_eltype),length(model))
    aux = zeros(@f(Base.promote_eltype),length(model) + 1 + B)
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
function ∂lnϕ∂n∂P(model::EoSModel, p, T, z=SA[1.], cache = ∂lnϕ_cache(model,p,T,z,Val{false}()); phase=:unknown, vol0=nothing)
    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
    RT = R̄*T
    V = volume(model, p, T, z, phase=phase, vol0=vol0)
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
function ∂lnϕ∂n∂P∂T(model::EoSModel, p, T, z=SA[1.],cache = ∂lnϕ_cache(model,p,T,z,Val{true}()); phase=:unknown, vol0=nothing)
    result,aux,lnϕ,∂lnϕ∂n,∂lnϕ∂P,∂P∂n,∂lnϕ∂T,hconfig = cache
    RT = R̄*T
    V = volume(model, p, T, z, phase=phase, vol0=vol0)
    n = sum(z)
    Z = p*V/RT/n

    ncomponents = length(z)
    aux[1] = V
    aux[2] = T
    aux[3:end] .= z
    F_res(model, V, T, z) = eos_res(model, V, T, z) / R̄ / T
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
