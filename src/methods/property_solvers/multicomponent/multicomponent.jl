
#original
#f(z) = eos(model,V,T,z)
#H(z) = ForwardDiff.hessian(f,z)/(R̄*T)
#L(z) = det(H(z))
#dL(z) = ForwardDiff.gradient(L,z)
#M(z) = [H(z)[1:end-1,:];transpose(dL(z))]
"""
    mixture_critical_constraint(model,V,T,z)

with `a(x)` the reduced `(A/RT)` helmholtz energy dependent on composition `xᵢ` for `i` ∈ `1:n`, returns `L` and `det(M)`, where `L` and `M` are defined as:
```
L := det(ℍ(a)) (ℍ = hessian)
M := ℍ(a) for rows ∈ 1:n-1
  := ∇L for row n
```
"""
function mixture_critical_constraint(model,V,T,z)
    f(x) = eos(model,V,T,x)/(Rgas(model)*T)
    H(x) = ForwardDiff.hessian(f,x) #∂A/∂zᵢ∂zⱼ == ∂A/∂zⱼ∂zᵢ
    L(x) = det(Symmetric(H(x)))
    dL(x) = ForwardDiff.gradient(L,x)
    HH = H(z)
    LL = det(HH)
    Mᵢ = @view(HH[end,:])
    Mᵢ .=  dL(z)
    MM = HH
    #M(x) = [HH[1:end-1,:];transpose(dL(x))]
    return LL , det(MM)
end

function μp_equality(model,v,T,w)
    np = length(v)
    nc = length(model)
    F = zeros(nc*(np - 1) + np - 1)
    return μp_equality(model, F, T, v, w)
end

function v_from_η(model::EoSModel, η, T, z)
    #we want a transformation such
    #v = lb -> η = -Inf
    #v = Inf, η = Inf
    #η = log(v - lb)
    #v = exp(η) + lb
   lb = lb_volume(model,T,z)
   V = exp(η) + lb/sum(z)
end

function v_from_η(model, model_r, η, T, z)
    if model_r == nothing
        return v_from_η(model, η, T, z)
    else
        return v_from_η(model_r, η, T, z)
    end
end

function η_from_v(model::EoSModel, V, T, z)
    lb =lb_volume(model,T,z)
    return log((V - lb)/sum(z))
end

function η_from_v(model::EoSModel,model_r, V, T, z)
    if model_r == nothing
        return η_from_v(model,V,T,z)
    else
        return η_from_v(model_r,V,T,z)
    end
end

struct TPspec{TT}
    T::TT
    p::TT
    pressure_specified::Bool
end

Tspec(T) = TPspec(T,zero(T)/zero(T),false)

function Pspec(p,T)
    _p,_T = promote(p,T)
    return TPspec(_T,_p,true)
end

function μp_equality(model::EoSModel, F, PT::TPspec, Base.@specialize(v), Base.@specialize(w))
    p,T = PT.p,PT.T
    R = Rgas(model)
    RTinv = 1/(R*T)
    w1,v1 = w[1],v[1]
    n_c = length(w1)
    n_p = length(v)

    μ1 = μj = similar(F,length(model))
    p1 = pressure(model,v1,T,w1)
    ∑w1 = sum(w1)
    VT_chemical_potential_res!(μ1,model,v1,T,w1)
    for j in 1:(n_p - 1)
        Fj = @inbounds viewn(F,n_c,j)
        for i in 1:n_c
            Fj[i] = μ1[i]
        end
    end

    p⁻¹ = 1/p_scale(model,w1)
    idx_p_start = n_c*(n_p - 1) + 1
    idx_p_end = n_c*(n_p - 1) + n_p - 1
    Fp = view(F,idx_p_start:idx_p_end)
    for j in 1:(n_p - 1)
        vj = @inbounds v[j+1]
        wj = @inbounds w[j+1]
        pj = pressure(model,vj,T,wj)
        #pᵣ_res_j = pressure_res(model,vj,T,wj)*RTinv
        #Δp = pᵣ_res_1 - pᵣ_res_j
        #Δp += pᵣ_ideal_1 - sum(wj)/vj
        Fp[j] = (p1 - pj)*p⁻¹
        VT_chemical_potential_res!(μ1,model,vj,T,wj)
        Fj = viewn(F,n_c,j)

        for i in 1:n_c
            μ1i = Fj[i]
            μji = μj[i]
            Δuᵣ = μ1i - μji
            Δu = Δuᵣ*RTinv + log(vj*w1[i]/(v1*wj[i]))
            Fj[i] = Δu
        end
    end

    if PT.pressure_specified
        #=
        p = -deos/dv = -deos_res/dt - d_eos_ideal/dv
        p = RT*(a_res)/dv - sum(z)*RT/v
        p/RT = -da_res/dv - 1/v
         =#
        #pᵣ_1 = pᵣ_res_1 + pᵣ_ideal_1
        F[idx_p_end + 1] = (p1 - p)*p⁻¹
    end

    return F
end

#non-condensable/non-volatile version
function μp_equality2(models::NTuple{2,M}, F, PT::TPspec, v, w, short_view) where M <: EoSModel
    p,T = PT.p,PT.T
    model_long, model_short = models
    v_long, v_short = v
    x_long, x_short = w
    n_short = length(x_short)
    n_long = length(x_long)
    μ_long = similar(F,n_long)
    μ_long = VT_chemical_potential_res!(μ_long,model_long,v_long,T,x_long)
    p_long = pressure(model_long,v_long,T,x_long)
    p_short = pressure(model_short,v_short,T,x_short)
    p⁻¹,RT⁻¹ = equilibria_scale(model_long,x_long)
    RTinv = 1/(Rgas(model_long)*T)
    μ_long_view = @view(μ_long[short_view])
    x_long_view = @view(x_long[short_view])
    for i in 1:n_short
        F[i] = μ_long_view[i]
    end
    μ_short = resize!(μ_long,n_short)
    μ_short = VT_chemical_potential_res!(μ_short,model_short,v_short,T,x_short)
    for i in 1:n_short
        μ_long_i = F[i]
        μ_short_i = μ_short[i]
        Δuᵣ = (μ_long_i - μ_short_i)
        Δμ = Δuᵣ*RTinv + log(v_short*x_long_view[i]/(v_long*x_short[i]))
        F[i] = Δμ*RT⁻¹
    end
    F[n_short+1] = (p_long-p_short)*p⁻¹
    if PT.pressure_specified
        F[n_short+2] = (p_long-p)*p⁻¹
    end
    return F
end

function μp_equality2(model::EoSModel,::Nothing, F, T, v, w, _view)
    return μp_equality(model,F,T,v,w)
end

function μp_equality2(model::EoSModel,model2::EoSModel, F, T, v, w, _view)
    return μp_equality2((model,model2),F,T,v,w,_view)
end

function wilson_k_values(model::EoSModel,p,T,crit = nothing)
    K = zeros(typeof(p+T+one(eltype(model))),length(model))
    return wilson_k_values!(K,model,p,T,crit)
end

function wilson_k_values!(K,model::EoSModel,p,T,crit = nothing)
    n = length(model)
    pure = split_model.(model)
    if crit === nothing
        crit = crit_pure.(pure)
    end
    for i ∈ 1:n
        pure_i = pure[i]
        Tc,pc,_ = crit[i]
        ps = first(saturation_pressure(pure_i,0.7*Tc))
        ω = -log10(ps/pc) - 1.0
        K[i] = exp(log(pc/p)+5.373*(1+ω)*(1-Tc/T))
    end
    return K
end

function bubbledew_check(vl,vv,zin,zout)
    (isapprox(vl,vv) && isapprox(zin,zout)) && return false
    !all(isfinite,zout) && return false
    !isfinite(vv) && return false
    !all(>=(0),zin) && return false
    !all(>=(0),zout) && return false
    return true
end

#generator for candidate fractions, given an initial composition, method by Pereira et al. (2010).
function initial_candidate_fractions(n)

    nc = length(n)
    x̂ = [zeros(nc) for i in 1:nc-1]
    x̄ = [zeros(nc) for i in 1:nc-1]

    for i ∈ 1:nc-1
        x̂i = x̂[i]
        x̄i = x̄[i]
        x̂i[i] = n[i]/2
        x̄i[i] = (1+n[i])/2
        for k ∈ 1:nc-1
            if k != i
                x̂i[k] = (1-x̂i[i])/(nc-1)
                x̄i[k] = (1-x̄i[i])/(nc-1)
            end
        end
        x̂i[nc] = 1 - sum(x̂i)
        x̄i[nc] = 1 - sum(x̄i)
    end
    x = vcat(x̂,x̄)
    return x
end

#when we have a candidate fraction, we generate points closer to that point.
function near_candidate_fractions(n,k = 0.5*minimum(n))
    nc = length(n)
    x = [zeros(nc) for i in 1:nc]
    for i in eachindex(x)
        xi = x[i]
        xi .= n
        xi[i] += k*n[i]
        xi ./= sum(xi)
    end
    return x
end

include("fugacity.jl")
include("rachford_rice.jl")
include("bubble_point.jl")
include("dew_point.jl")
include("azeotrope_point.jl")
include("LLE_point.jl")
include("VLLE.jl")
include("crit_mix.jl")
include("UCEP.jl")
include("UCST_mix.jl")
include("tp_flash.jl")
include("krichevskii_parameter.jl")
include("solids/sle_solubility.jl")
include("solids/slle_solubility.jl")
include("solids/eutectic_point.jl")

export bubble_pressure_fug, bubble_temperature_fug, dew_temperature_fug, dew_pressure_fug
export bubble_pressure,    dew_pressure,    LLE_pressure,    azeotrope_pressure, VLLE_pressure
export bubble_temperature, dew_temperature, LLE_temperature, azeotrope_temperature, VLLE_temperature
export crit_mix, UCEP_mix, UCST_mix
export krichevskii_parameter
export sle_solubility, eutectic_point, slle_solubility
