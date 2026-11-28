function __psi_res(model, T, rhovec)
    nc = length(model)
    _rhovec = @view rhovec[1:nc]
    return psi_res(model, T, _rhovec)
end

function __psi_res(model, data)
    nc = length(model)
    rhovec = @view data[1:nc]
    T = data[nc + 1]
    return psi_res(model, T, rhovec)
end

StaticForwardDiffTags.inner_function(::typeof(__psi_res)) = psi_res

#returns μ_res/RT and p_res/RT from isochoric derivatives of Ψ = A_res/(RT·V), V = 1.
function __μpr(model, T, ρ, cache)
    nc = length(model)

    if is_idealmodel(model)
        TT = Base.promote_eltype(model, T, ρ)
        if cache isa ∂lnϕCache
            μr = cache.lnϕ
            fill!(μr, 0)
        elseif cache isa AbstractVector
            μr = cache
            fill!(μr, 0)
        else
            nz = length(ρ)
            μr = zeros(TT, nz)
        end
        return μr, zero(TT)*one(eltype(μr))
    end

    #deferred forwarddiff tag for faster AD in some cases
    Ψ = @deferred_isochoric_rhovec(__psi_res, psi_res)
    if cache isa ∂lnϕCache
        aux = cache.aux
        ρᵢ = @view aux[1:nc]
        ρᵢ .= ρ
        gconfig = Solvers._GradientConfig(cache.hconfig)
        _result = ForwardDiff.gradient!(cache.result, Ψ, aux, gconfig, Val{false}())
    else
        if cache isa AbstractVector
            result = DiffResults.MutableDiffResult(zero(eltype(cache)), (cache,))
        else
            result = DiffResults.GradientResult(ρ)
        end
        _result = ForwardDiff.gradient!(result, Ψ, ρ)
        ρᵢ = ρ
    end

    μ_res = @view DiffResults.gradient(_result)[1:nc]

    p_res = dot(μ_res, ρᵢ) - DiffResults.value(_result)
    return μ_res, p_res
end

# Function to compute fugacity coefficient
function lnϕ(model::EoSModel, p, T, z=SA[1.0], cache=nothing; phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    ∑z = sum(z)
    RT = Rgas(model)*T
    logZ = log(p*vol/RT/∑z)
    nc = length(z)

    if cache isa Vector
        return lnϕ!(cache, model, p, T, z; vol)
    elseif cache isa ∂lnϕCache
        lnϕ_out = cache.lnϕ
        if nc == 1
            lnϕ_out[1] = VT_∑zlogϕ(model, vol/∑z, T, SA[1.0], p)
        else
            rhovec = @view cache.aux[1:nc]
            rhovec .= z ./ vol
            μ_res_RT, _ = __μpr(model, T, rhovec, cache)
            lnϕ_out .= μ_res_RT .- logZ
        end
    else
        μ_res = VT_chemical_potential_res(model, vol, T, z)
        if ismutable(μ_res)
            lnϕ_out = μ_res
            lnϕ_out .= μ_res ./ RT .- logZ
        else
            lnϕ_out = μ_res/RT .- logZ
        end
    end
    return lnϕ_out, vol
end

function lnϕ(model::IdealModel, p, T, z=SA[1.0], cache=nothing; phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    lnϕ = FillArrays.Zeros(length(z))
    return lnϕ, vol
end

function lnϕ!(cache::∂lnϕCache, model::EoSModel, p, T, z=SA[1.0]; phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    is_idealmodel(model) && return (fill!(cache[3], 0.0), vol)

    return lnϕ(model, p, T, z, cache; vol)
end

function lnϕ!(cache::Nothing, model::EoSModel, p, T, z=SA[1.0]; phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    return lnϕ(model, p, T, z; vol)
end

function lnϕ!(lnϕ_vec::AbstractVector, model::EoSModel, p, T, z=SA[1.0], cache=nothing; phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    if is_idealmodel(model)
        lnϕ_vec .= 0
        return lnϕ_vec, vol
    end

    if cache !== nothing
        lnϕx, _, = lnϕ(model, p, T, z, cache; vol)
        lnϕ_vec .= lnϕx
    else
        RT = Rgas(model)*T
        μ_res = VT_chemical_potential_res!(lnϕ_vec, model, vol, T, z)
        Z = p*vol/RT/sum(z)
        lnϕ_vec .= μ_res ./ RT .- log(Z)
    end
    return lnϕ_vec, vol
end

function ∑zlogϕ(model::EoSModel, p, T, z=SA[1.0], cache=nothing; phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    return VT_∑zlogϕ(model, vol, T, z, p), vol
end

#ideal models: ψ_res = 0, so μ_res, p_res, H and ∂(p_res/RT)/∂T are zero. nH is the size of H.
function __μpr_hess_ideal(nc, T, rhovec, cache, nH)
    TT = Base.promote_eltype(model, T, rhovec)
    nc = length(model)
    if cache === nothing
        μ_res = zeros(TT, nc)
        H = zeros(TT, (nH, nH))
    else
        μ_res = cache.lnϕ
        μ_res .= 0
        result = cache.result
        H = result.derivs[2]
        H .= 0
    end

    return μ_res, zero(TT), H, zero(TT)
end

#μ_res/RT, p_res/RT and the Hessian of Ψ over the nc densities.
function __μpr_hess(model, T, rhovec, cache)
    nc = length(rhovec)

    if is_idealmodel(model)
        Tx = T*one(eltype(model))
        return __μpr_hess_ideal(nc, Tx, rhovec, cache, nc)
    end

    Ψ = @deferred_isochoric_rhovec(__psi_res, psi_res)

    if cache === nothing
        result = DiffResults.HessianResult(rhovec)
        _result = ForwardDiff.hessian!(result, Ψ, rhovec)
    else
        aux = cache.aux
        rhovec_aux = @view aux[1:nc]
        rhovec_aux .= rhovec
        _result = ForwardDiff.hessian!(cache.result, Ψ, aux, cache.hconfig, Val{false}())
    end

    μ_res = @view DiffResults.gradient(_result)[1:nc]
    H = @view DiffResults.hessian(_result)[1:nc, 1:nc]
    p_res = dot(μ_res, rhovec) - DiffResults.value(_result)
    return μ_res, p_res, H, zero(p_res)
end

#μ_res/RT, p_res/RT, the Hessian of Ψ, and ∂(p_res/RT)/∂T; over the nc densities.
function __μprT_hess(model, T, rhovec, cache)
    nc = length(rhovec)

    if is_idealmodel(model)
        Tx = T*one(eltype(model))
        return __μpr_hess_ideal(nc, Tx, rhovec, cache, nc + 1)
    end

    Ψ = WithContext(Base.Fix1(__psi_res, model), ∂Tag{psi_res})

    if cache === nothing
        aux = vcat(rhovec, T)
        result = DiffResults.HessianResult(aux)
        _result = ForwardDiff.hessian!(result, Ψ, aux)
    else
        aux = cache.aux
        rhovec_aux = @view aux[1:nc]
        rhovec_aux .= rhovec
        aux[nc + 1] = T
        _result = ForwardDiff.hessian!(cache.result, Ψ, aux, cache.hconfig, Val{false}())
    end

    ∇Ψ = DiffResults.gradient(_result)
    μ_res = @view ∇Ψ[1:nc]
    H = @view DiffResults.hessian(_result)[1:(nc + 1), 1:(nc + 1)]
    p_res = dot(μ_res, rhovec) - DiffResults.value(_result)
    dpdT = -∇Ψ[nc + 1]
    for i in 1:nc
        dpdT += rhovec[i]*H[i, nc + 1]
    end
    return μ_res, p_res, H, dpdT
end

#fills lnϕ and a (in the cache, if available). returns p/RT, a·ρ, the Hessian view, ∂(p/RT)/∂T|ρ.
function __lnϕ_core!(model, T, rhovec, cache, mup_hess_fn::F) where F
    μ, p_res, H, dptotdT = mup_hess_fn(model, T, rhovec, cache)
    ρ = sum(rhovec)
    lnZ = log1p(p_res/ρ)
    if cache === nothing
        a = similar(μ)
        if ismutable(μ)
            μ .= μ .- lnZ
            lnϕ = μ
        else
            lnϕ = μ .- lnZ
        end
    else
        a = cache.scratch
        lnϕ = cache.lnϕ
        lnϕ .= μ .- lnZ
    end

    #a = ∂(p/RT)/∂ρ = 1 + Hρ. returns a·ρ.
    aρ = zero(eltype(a))
    for i in eachindex(rhovec)
        ai = one(eltype(a))
        for k in eachindex(rhovec)
            ai += H[i, k]*rhovec[k]
        end
        a[i] = ai
        aρ += ai*rhovec[i]
    end

    return p_res, aρ, H, lnϕ, a, dptotdT
end

_dlnϕdT!(cache::∂lnϕCache, a, H, T, p_over_RT, aρ, dptotdT) = _dlnϕdT!(cache.∂lnϕ∂T, a, H, T, p_over_RT, aρ, dptotdT)
_dlnϕdT!(::Nothing, a, H, T, p_over_RT, aρ, dptotdT) = _dlnϕdT!(similar(a), a, H, T, p_over_RT, aρ, dptotdT)

function _dlnϕdT!(∂lnϕ∂T, a, H, T, p_over_RT, aρ, dptotdT)
    nc = length(a)
    s = -(p_over_RT/T + dptotdT)/aρ
    for i in 1:nc
        ∂lnϕ∂T[i] = H[i, nc + 1] + 1/T + s*a[i]
    end
    return ∂lnϕ∂T
end

_dlnϕdnP!(cache::∂lnϕCache, a, H, RT, p_over_RT, aρ, vol, n) = _dlnϕdnP!(cache.∂lnϕ∂n, cache.∂lnϕ∂P, a, H, RT, p_over_RT, aρ, vol, n)

function _dlnϕdnP!(::Nothing, a, H, RT, p_over_RT, aρ, vol, n)
    nc = length(a)
    ∂lnϕ∂n = similar(a, (nc, nc))
    ∂lnϕ∂P = similar(a)
    return _dlnϕdnP!(∂lnϕ∂n, ∂lnϕ∂P, a, H, RT, p_over_RT, aρ, vol, n)
end

function _dlnϕdnP!(∂lnϕ∂n, ∂lnϕ∂P, a, H, RT, p_over_RT, aρ, vol, n)
    ninv = 1/n
    nc = length(a)
    for i in 1:nc
        ai_aρ = a[i]/aρ
        ∂lnϕ∂P[i] = (ai_aρ - 1/p_over_RT)/RT
        for j in 1:nc
            ∂lnϕ∂n[i, j] = (H[i, j] - a[j]*ai_aρ)/vol + ninv
        end
    end
    return ∂lnϕ∂n, ∂lnϕ∂P
end

# Function to compute fugacity coefficient and its pressure and composition derivatives
function ∂lnϕ∂n∂P(model::EoSModel, p, T, z=SA[1.0], cache=∂lnϕ_cache(model, p, T, z, Val{false}()); phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    if is_idealmodel(model)
        if cache isa ∂lnϕCache
            fill!(cache.lnϕ, 0)
            fill!(cache.∂lnϕ∂n, 0)
            fill!(cache.∂lnϕ∂P, 0)
            return cache.lnϕ, cache.∂lnϕ∂n, cache.∂lnϕ∂P, vol
        else
            TT = Base.promote_eltype(model, p, T, z)
            nz = length(z)
            return zeros(TT, nz), zeros(TT, nz, nz), zeros(TT, nz), vol
        end
    end

    nc = length(cache.scratch)
    if cache === nothing
        rhovec = z ./ vol
    else
        aux = cache.aux
        rhovec = @view aux[1:nc]
        rhovec .= z ./ vol
    end

    n = sum(z)
    RT = Rgas(model)*T
    p_over_RT = p/RT
    _, aρ, H, lnϕ, a, _ = __lnϕ_core!(model, T, rhovec, cache, __μpr_hess)
    ∂lnϕ∂n, ∂lnϕ∂P = _dlnϕdnP!(cache, a, H, RT, p_over_RT, aρ, vol, n)
    return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, vol
end

# Function to compute fugacity coefficient pressure derivatives
function ∂lnϕ∂P(model::EoSModel, p, T, z=SA[1.0], cache=∂lnϕ_cache(model, p, T, z, Val{false}()); phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    if is_idealmodel(model)
        if cache isa ∂lnϕCache
            return fill!(cache.∂lnϕ∂P, 0), vol
        else
            TT = Base.promote_eltype(model, p, T, z)
            nz = length(z)
            return zeros(TT, nz), vol
        end
    end

    nc = length(cache.scratch)
    RT = Rgas(model)*T
    pᵣ(x) = pressure_res(model, 1.0, T, @view(x[1:nc]))/RT
    ρt = sum(z)/vol

    if cache === nothing
        rhovec = zeros(Base.promote_eltype(model, p, T, z), nc)
        rhovec .= z ./ vol
        _result = ForwardDiff.gradient!(DiffResults.GradientResult(rhovec), pᵣ, rhovec)
        ∇pᵣ = @view DiffResults.gradient(_result)[1:nc]
        aρ = ρt + dot(rhovec, ∇pᵣ)
        ∂lnϕ∂P = rhovec
    else
        aux = cache.aux
        rhovec = @view aux[1:nc]
        rhovec .= z ./ vol
        gconfig = Solvers._GradientConfig(cache.hconfig)
        _result = ForwardDiff.gradient!(cache.result, pᵣ, aux, gconfig, Val{false}())
        ∇pᵣ = @view DiffResults.gradient(_result)[1:nc]
        aρ = ρt + dot(rhovec, ∇pᵣ)
        ∂lnϕ∂P = cache.∂lnϕ∂P
    end

    p_over_RT = p/RT
    for i in 1:nc
        ∂lnϕ∂P[i] = ((1 + ∇pᵣ[i])/aρ - 1/p_over_RT)/RT
    end

    return ∂lnϕ∂P, vol
end

# Function to compute fugacity coefficient temperature derivatives
function ∂lnϕ∂T(model::EoSModel, p, T, z=SA[1.0], cache=∂lnϕ_cache(model, p, T, z, Val{false}()); phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    if is_idealmodel(model)
        if cache isa ∂lnϕCache
            return fill!(cache.∂lnϕ∂P, 0), vol
        else
            TT = Base.promote_eltype(model, p, T, z)
            nz = length(z)
            return zeros(TT, nz), vol
        end
    end

    if cache === nothing
        if is_idealmodel(model)
        end
        rhovec = z ./ vol
    else
        nz = length(z)
        rhovec = @view cache.aux[1:nz]
        rhovec .= z ./ vol
    end
    nc = length(z)

    s = -VT_isobaric_expansivity(model, vol, T, z)

    #T-derivative of ρₜ·a_res along the isobar, differentiated w.r.t. ρ
    function Φ(x)
        ρ = @view x[1:nc]
        ρt = sum(ρ)
        f(ε) = ρt*a_res(model, 1/(1 + s*ε), T + ε, ρ)
        return Solvers.derivative(f, zero(T))
    end

    if cache === nothing
        ∂lnϕ∂T = ForwardDiff.gradient(Φ, rhovec)
        ∂lnϕ∂T .+= 1/T + s
    else
        gconfig = Solvers._GradientConfig(cache.hconfig)
        _result = ForwardDiff.gradient!(cache.result, Φ, cache.aux, gconfig, Val{false}())
        ∂lnϕ∂T = cache.∂lnϕ∂T
        ∂lnϕ∂T .= @view(DiffResults.gradient(_result)[1:nc])
        ∂lnϕ∂T .+= 1/T + s
    end
    return ∂lnϕ∂T, vol
end

# Function to compute fugacity coefficient and its temperature, pressure and composition derivatives
function ∂lnϕ∂n∂P∂T(model::EoSModel, p, T, z=SA[1.0], cache=∂lnϕ_cache(model, p, T, z, Val{true}()); phase=:unknown, vol0=nothing, threaded=true, vol=volume(model, p, T, z; phase, vol0, threaded))
    if is_idealmodel(model)
        nz = length(z)
        if cache isa ∂lnϕCache
            fill!(cache.lnϕ, 0)
            fill!(cache.∂lnϕ∂n, 0)
            fill!(cache.∂lnϕ∂P, 0)
            if length(cache.aux) == nz
                ∂lnϕ∂T0 = cache.scratch
            else
                ∂lnϕ∂T0 = cache.∂lnϕ∂T
            end
            fill!(∂lnϕ∂T0, 0)
            return cache.lnϕ, cache.∂lnϕ∂n, cache.∂lnϕ∂P, ∂lnϕ∂T0, vol
        else
            TT = Base.promote_eltype(model, p, T, z)
            return zeros(TT, nz), zeros(TT, nz, nz), zeros(TT, nz), zeros(TT, nz), vol
        end
    end

    nc = length(model)
    if cache === nothing
        rhovec = z ./ vol
    else
        aux = cache.aux
        rhovec = @view aux[1:nc]
        rhovec .= z ./ vol
    end

    if cache isa ∂lnϕCache && nc == length(cache.aux)
        #not enough cache for one pass, so use two
        lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, vol = ∂lnϕ∂n∂P(model, p, T, z, cache; vol)
        ∂lnϕ∂T_out = cache.scratch
        ∂lnϕ∂T_out .= lnϕ                       #save lnϕ: cache.∂lnϕ∂T aliases cache.lnϕ here
        ∂lnϕ∂T(model, p, T, z, cache; vol)           #overwrites cache.lnϕ with ∂lnϕ/∂T
        for i in 1:nc                           #swap, so lnϕ holds lnϕ and scratch holds ∂lnϕ/∂T
            lnϕ[i], ∂lnϕ∂T_out[i] = ∂lnϕ∂T_out[i], lnϕ[i]
        end
        return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, ∂lnϕ∂T_out, vol
    end

    n = sum(z)
    RT = Rgas(model)*T
    p_over_RT = p/RT
    _, aρ, H, lnϕ, a, dptotdT = __lnϕ_core!(model, T, rhovec, cache, __μprT_hess)
    ∂lnϕ∂n, ∂lnϕ∂P = _dlnϕdnP!(cache, a, H, RT, p_over_RT, aρ, vol, n)
    ∂lnϕ∂T = _dlnϕdT!(cache, a, H, T, p_over_RT, aρ, dptotdT)
    return lnϕ, ∂lnϕ∂n, ∂lnϕ∂P, ∂lnϕ∂T, vol
end

#functions used to overload gamma-phi models
function modified_lnϕ(model, p, T, z, cache; phase=:unknown, vol0=nothing)
    lnϕz, vz = lnϕ(model, p, T, z, cache; phase, vol0)
    if isnan(vz)
        lnϕz, vz = lnϕ(model, p, T, z, cache; phase)
    end
    return lnϕz, vz
end

modified_gibbs(model, p, T, w) = modified_gibbs(model, p, T, w, :unknown, oftype(zero(Base.promote_eltype(model, p, T, w)), NaN))
modified_gibbs(model, p, T, w, phase) = modified_gibbs(model, p, T, w, phase, oftype(zero(Base.promote_eltype(model, p, T, w)), NaN))

function modified_gibbs(model, p, T, w, phase, vol)
    if isnan(vol) || isnothing(vol)
        volw = volume(model, p, T, w, phase=phase)
    else
        volw = vol
    end
    RT = Rgas(model)*T
    g = PT_property(model, p, T, w, phase, volw, VT_gibbs_energy)
    return g/RT, volw
end

_dlnϕdn!(cache::∂lnϕCache, a, H, RT, p_over_RT, aρ, vol, n) = _dlnϕdn!(cache.∂lnϕ∂n, a, H, RT, p_over_RT, aρ, vol, n)

function _dlnϕdn!(::Nothing, a, H, RT, p_over_RT, aρ, vol, n)
    nc = length(a)
    ∂lnϕ∂n = similar(a, (nc, nc))
    return _dlnϕdn!(∂lnϕ∂n, a, H, RT, p_over_RT, aρ, vol, n)
end

function _dlnϕdn!(∂lnϕ∂n, a, H, RT, p_over_RT, aρ, vol, n)
    ninv = 1/n
    nc = length(a)
    for i in 1:nc
        ai_aρ = a[i]/aρ
        for j in 1:nc
            ∂lnϕ∂n[i, j] = (H[i, j] - a[j]*ai_aρ)/vol + ninv
        end
    end
    return ∂lnϕ∂n
end

function modified_∂lnϕ∂n(model, p, T, z, cache; phase=:unknown, vol0=nothing)
    vol = volume(model, p, T, z; phase, vol0, threaded=false)

    if is_idealmodel(model)
        nz = length(z)
        if cache isa ∂lnϕCache
            fill!(cache.lnϕ, 0)
            fill!(cache.∂lnϕ∂n, 0)
            return cache.lnϕ, cache.∂lnϕ∂n, vol
        else
            TT = Base.promote_eltype(model, p, T, z)
            return zeros(TT, nz), zeros(TT, nz, nz), vol
        end
    end

    nc = length(model)
    if cache === nothing
        rhovec = z ./ vol
    else
        aux = cache.aux
        rhovec = @view aux[1:nc]
        rhovec .= z ./ vol
    end

    n = sum(z)
    RT = Rgas(model)*T
    p_over_RT = p/RT
    _, aρ, H, lnϕ, a, _ = __lnϕ_core!(model, T, rhovec, cache, __μpr_hess)
    ∂lnϕ∂n = _dlnϕdn!(cache, a, H, RT, p_over_RT, aρ, vol, n)
    return lnϕ, ∂lnϕ∂n, vol
end

function lnf(model, V, T, z, cache=nothing)
    RT = Rgas(model)*T
    n = sum(z)
    logZp = log(V/RT/n)
    nc = length(z)
    if cache isa ∂lnϕCache
        lnf = cache.lnϕ
        if nc == 1
            lnf1, p = VT_lnf_pure(model, V, T)
            lnf[1] = lnf1
        else
            rhovec = @view cache.aux[1:nc]
            rhovec .= z ./ V
            μ_res_RT, p_res_RT = __μpr(model, T, rhovec, cache)
            p = RT*(p_res_RT + n/V)
            lnf .= μ_res_RT .- logZp
        end
    else
        rhovec = z ./ V
        μ_res_RT, p_res_RT = __μpr(model, T, rhovec, cache)
        p = RT*(p_res_RT + n/V)
        if ismutable(μ_res_RT)
            lnf = μ_res_RT
            lnf .= μ_res_RT .- logZp
        else
            lnf = μ_res_RT .- logZp
        end
    end
    return lnf, p
end

function VT_lnf_pure(model, V, T)
    RT = Rgas(model)*T
    F, dFdV = f∂fdV_res(model, V, T, SA[1.0])
    p_res = -dFdV
    μ_res = F + p_res*V
    Zp = V/RT
    lnf = μ_res/RT - log(Zp)
    p = p_res + RT/V
    return lnf, p
end
