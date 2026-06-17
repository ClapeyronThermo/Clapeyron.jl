
#original
#f(z) = eos(model,V,T,z)
#H(z) = ForwardDiff.hessian(f,z)/(R̄*T)
#L(z) = det(H(z))
#dL(z) = ForwardDiff.gradient(L,z)
#M(z) = [H(z)[1:end-1,:];transpose(dL(z))]
"""
    mixture_critical_constraint(model,V,T,z)

with `a(x)` the reduced `(A/RT)` Helmholtz energy dependent on composition `xᵢ` for `i` ∈ `1:n`, returns `L` and `det(M)`, where `L` and `M` are defined as:
```
L := det(ℍ(a)) (ℍ = hessian)
M := ℍ(a) for rows ∈ 1:n-1
  := ∇L for row n
```
"""
function mixture_critical_constraint(model,V,T,z)
    f(x) = sum(x)*(a_res(model,V,T,x) + a_ideal(BasicIdeal(),V,T,x))
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

η_from_v(model::EoSModel,::Nothing, V, T, z) = η_from_v(model,V,T,z)
η_from_v(model::EoSModel,model_r::EoSModel, V, T, z) = η_from_v(model_r,V,T,z)


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
    VT_chemical_potential_res!(μ1,model,v1,T,w1)
    log_v1 = log(v1)
    @inbounds for j in 1:(n_p - 1)
        Fj = viewn(F,n_c,j)
        for i in 1:n_c
            Fj[i] = μ1[i]
        end
    end

    p⁻¹ = 1/p_scale(model,w1)
    idx_p_start = n_c*(n_p - 1) + 1
    idx_p_end = n_c*(n_p - 1) + n_p - 1
    Fp = view(F,idx_p_start:idx_p_end)
    @inbounds for j in 1:(n_p - 1)
        vj = v[j+1]
        wj = w[j+1]
        pj = pressure(model,vj,T,wj)
        #pᵣ_res_j = pressure_res(model,vj,T,wj)*RTinv
        #Δp = pᵣ_res_1 - pᵣ_res_j
        #Δp += pᵣ_ideal_1 - sum(wj)/vj
        Fp[j] = (p1 - pj)*p⁻¹
        VT_chemical_potential_res!(μ1,model,vj,T,wj)
        Fj = viewn(F,n_c,j)

        log_v_common = log(vj) - log_v1
        for i in 1:n_c
            μ1i = Fj[i]
            μji = μj[i]
            Δuᵣ = μ1i - μji
            Fj[i] = Δuᵣ*RTinv + log_v_common + log(w1[i]) - log(wj[i])
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
    @inbounds for i in 1:n_short
        F[i] = μ_long_view[i]
    end
    μ_short = resize!(μ_long,n_short)
    log_v_common = log(v_short/v_long)
    if n_short == 1
        ∑n_short = sum(x_short)
        p_res = p_short - ∑n_short*Rgas(model_short)*T/v_short
        μ_short_1 = (eos_res(model_short,v_short,T,x_short) + p_res*v_short)/∑n_short
        μ_long_1 = F[1]
        Δuᵣ = μ_long_1 - μ_short_1
        Δμ = Δuᵣ*RTinv + log_v_common + (log(x_long_view[1]) - log(x_short[1]))
        F[1] = Δμ*RT⁻¹
    else
        VT_chemical_potential_res!(μ_short,model_short,v_short,T,x_short)
        @inbounds for i in 1:n_short
            μ_long_i  = F[i]
            μ_short_i = μ_short[i]
            Δuᵣ = μ_long_i - μ_short_i
            Δμ = Δuᵣ*RTinv + log_v_common + (log(x_long_view[i]) - log(x_short[i]))
            F[i] = Δμ*RT⁻¹
        end
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

wilson_k_values!(K,model::EoSModel,p,T) = wilson_k_values!(K,model,p,T,nothing)

function wilson_k_values!(K,model::EoSModel,p,T,crit)
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
        K[i] = exp(log(pc/p)+5.3726985503194395*(1+ω)*(1-Tc/T)) #5.37 = log(10)*7/3
    end
    return K
end

function bubbledew_check(model,p,T,vw,vz,w,z)
    vmin,vmax = minmax(vw,vz)
    dv = (vmax - vmin)/vmax
    dz = z_norm(z,w)
    dz < 1e-5 && dv < 1e-3 && return false
    dv < 1e-5 && dz < 1e-3 && return false
    !all(isfinite,w) && return false
    !isfinite(vw) && return false
    !all(>=(0),w) && return false
    !all(>=(0),z) && return false
    if has_a_res(model) && !(any(iszero,w)) #the second check is to exclude nonvolatiles/noncondensables. TODO: find a better way to do this.
        #all normal checks are ok, now we check if the origin phase result is really the most stable one
        gz = VT_gibbs_free_energy(model,vz,T,z)
        gz_p = gibbs_free_energy(model,p,T,z)
        gz_w = VT_gibbs_free_energy(model,vw,T,w)
        dg = (gz-gz_p)/gz
        if gz_p < gz && abs(dg) > 0.001
            return false
        end
    end
    return true
end

function initial_candidate_fractions(n::AbstractVector{TT}) where TT
    x = Vector{TT}[]
    nc = length(n)

    for i in 1:nc-1
        push!(x,z_pereira!(similar(n),n,i,true))
    end
    for i in 1:nc-1
        push!(x,z_pereira!(similar(n),n,i,false))
    end
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

function bubbledew_pressure_ad_v(result,tup,λtup,_bubble)
    f(x,tups) = begin
        model,T,z = tups
        vl = x[1]
        vv = x[2]
        w = @view x[3:end]
        if _bubble
            _x,_y = z,w
        else
            _x,_y = w,z
        end
        lnfl,pl = lnf(model,vl,T,_x)
        lnfv,pv = lnf(model,vv,T,_y)
        
        F1 = pl - pv
        F2 = sum(w) - 1.0 # can exclude this restriction, but would then need additional logic to parse w (excluding one component)
        F3 = lnfl - lnfv + log.(_x) - log.(_y)
        res = vcat(F1,F2,F3) # can probably be efficient with preallocation and @view but requires the common Dual type between tups and x, otherwise __gradients_for_root_finders will have the incorrect Dual type
        return res
    end
    λx = vcat(result[2],result[3],result[4])
    ∂x = __gradients_for_root_finders(λx,tup,λtup,f)
    ∂vl,∂vv = ∂x[1],∂x[2]
    ∂model,∂T,∂z = tup
    ∂p = _bubble ? pressure(∂model,∂vl,∂T,∂z) : pressure(∂model,∂vv,∂T,∂z)
    ∂w = ∂x[3:end]
    return ∂p,∂vl,∂vv,∂w
end

function bubbledew_pressure_ad_p(result,tup,λtup,_bubble,lle = false)
    vl0,vv0 = result[2],result[3]
    f(x,tups) = begin
        model,T,z = tups
        p = x[1]
        w = @view x[2:end]
        if _bubble
            _x,_y = z,w
        else
            _x,_y = w,z
        end

        phasey = lle ? :liquid : :vapour

        lnϕl,_ = modified_lnϕ(model,p,T,_x,nothing,phase = :liquid,vol0 = primalval(vl0))
        lnϕv,_ = modified_lnϕ(model,p,T,_y,nothing,phase = phasey, vol0 = primalval(vv0))
        
        F1 = sum(w) - 1.0 # can exclude this restriction, but would then need additional logic to parse w (excluding one component)
        F2 = lnϕl - lnϕv + log.(_x) - log.(_y)
        res = vcat(F1,F2) # can probably be efficient with preallocation and @view but requires the common Dual type between tups and x, otherwise __gradients_for_root_finders will have the incorrect Dual type
        return res
    end
    λp = result[1]
    λw = result[4]
    λmodel,λT,λz = λtup
    ∂model,∂T,∂z = tup
    λx = vcat(λp,λw)
    ∂x = __gradients_for_root_finders(λx,tup,λtup,f)
    ∂p = ∂x[1]
    ∂w = ∂x[2:end]
    #∂vl = volume_ad(result[2],(∂model,∂p,∂T,∂z),(λmodel,λp,λT,λz))
    if _bubble
        _∂x,_∂y = ∂z,∂w
    else
        _∂x,_∂y = ∂w,∂z
    end
    #∂vv = volume_ad(result[3],(∂model,∂p,∂T,∂w),(λmodel,λp,λT,λy))
    ∂vl = volume(∂model,∂p,∂T,_∂x,phase = :liquid,vol0 = primalval(vl0))
    phasey = lle ? :liquid : :vapour
    ∂vv = volume(∂model,∂p,∂T,_∂y,phase = phasey,vol0 = primalval(vv0))
    return ∂p,∂vl,∂vv,∂w
end

function bubbledew_temperature_ad_v(result,tup,λtup,_bubble)
    f(x,tups) = begin
        model,p,z = tups
        T = x[1]
        vl = x[2]
        vv = x[3]
        w = @view x[4:end]
        if _bubble
            _x,_y = z,w
        else
            _x,_y = w,z
        end
        lnfl,pl = lnf(model,vl,T,_x)
        lnfv,pv = lnf(model,vv,T,_y)
        F1 = pl - p
        F2 = pv - p
        F3 = sum(w) - 1.0 # can exclude this restriction, but would then need additional logic to parse w (excluding one component)
        F4 = lnfl - lnfv + log.(_x) - log.(_y)
        vcat(F1,F2,F3,F4) # can probably be efficient with preallocation and @view but requires the common Dual type between tups and x, otherwise __gradients_for_root_finders will have the incorrect Dual type
    end
    λx = vcat(result[1],result[2],result[3],result[4])
    ∂x = __gradients_for_root_finders(λx,tup,λtup,f)
    ∂T,∂vl,∂vv = ∂x[1:3]
    ∂w = ∂x[4:end]
    return ∂T,∂vl,∂vv,∂w
end

function bubbledew_temperature_ad_p(result,tup,λtup,_bubble,lle)
    vl0,vv0 = result[2],result[3]
    f(x,tups) = begin
        model,p,z = tups
        T = x[1]
        w = @view x[2:end]
        if _bubble
            _x,_y = z,w
        else
            _x,_y = w,z
        end

        phasey = lle ? :liquid : :vapour

        lnϕl,_ = modified_lnϕ(model,p,T,_x,nothing,phase = :liquid,vol0 = primalval(vl0))
        lnϕv,_ = modified_lnϕ(model,p,T,_y,nothing,phase = phasey, vol0 = primalval(vv0))
        F1 = sum(w) - 1.0 # can exclude this restriction, but would then need additional logic to parse w (excluding one component)
        F2 = lnϕl - lnϕv + log.(_x) - log.(_y)
        vcat(F1,F2) # can probably be efficient with preallocation and @view but requires the common Dual type between tups and x, otherwise __gradients_for_root_finders will have the incorrect Dual type
    end

    λT = result[1]
    λy = result[4]
    λmodel,λp,λz = λtup
    ∂model,∂p,∂z = tup
    λx = vcat(λT,λy)
    ∂x = __gradients_for_root_finders(λx,tup,λtup,f)
    ∂T = ∂x[1]
    ∂w = ∂x[2:end]
    ∂vl = volume_ad(vl0,(∂model,∂p,∂T,∂z),(λmodel,λp,λT,λz))
    ∂vv = volume_ad(vv0,(∂model,∂p,∂T,∂w),(λmodel,λp,λT,λy))
    return ∂T,∂vl,∂vv,∂w
end
 
bubble_temperature_ad(result,tup,λtup) = bubbledew_temperature_ad_v(result,tup,λtup,true)
dew_temperature_ad(result,tup,λtup) = bubbledew_temperature_ad_v(result,tup,λtup,false)
bubble_pressure_ad(result,tup,λtup) = bubbledew_pressure_ad_v(result,tup,λtup,true)
dew_pressure_ad(result,tup,λtup) = bubbledew_pressure_ad_v(result,tup,λtup,false)

function zero_non_equilibria!(w,in_equilibria)
    for i in eachindex(w)
        in_equilibria[i] || (w[i] = 0)
    end
    return w
end

function comps_in_equilibria(components,::Nothing)::Vector{Bool}
    return fill(true,length(components))
end

function comps_in_equilibria(components,not_in_w)::Vector{Bool}
    res = fill(true,length(components))
    for i in eachindex(res)
        res[i] = !in(components[i],not_in_w)
    end
    return res
end

include("fugacity.jl")
include("rachford_rice.jl")
include("bubble_point.jl")
include("dew_point.jl")
include("LLE_point.jl")
include("azeotrope_point.jl")
include("VLLE.jl")
include("crit_mix.jl")
include("UCEP.jl")
include("UCST_mix.jl")
include("flash.jl")

#include("bubbledew/ChemPotQX.jl")
#include("bubbledew/FugQX.jl")
#include("bubbledew/ActivityQT.jl")

include("krichevskii_parameter.jl")
include("solids/sle_solubility.jl")
include("solids/slle_solubility.jl")
include("solids/eutectic_point.jl")

export bubble_pressure,    dew_pressure,    LLE_pressure,    azeotrope_pressure, VLLE_pressure
export bubble_temperature, dew_temperature, LLE_temperature, azeotrope_temperature, VLLE_temperature
export crit_mix, UCEP_mix, UCST_pressure, UCST_temperature, UCST_mix, mechanical_critical_point
export krichevskii_parameter
export sle_solubility, sle_solubility_T, eutectic_point, slle_solubility