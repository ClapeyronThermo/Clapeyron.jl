
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

function μp_equality(model,v_l,v_v,T,x,y)
    F = zeros(length(model)+1)
    ps = p_scale(model,x)
    return μp_equality(model, F, T, v_l, v_v, x, y,ps)
end

function μp_equality(model::EoSModel, F, T, v_l, v_v, x, y,ps,Ts = T)
    n_c = length(x)
    p_l = pressure(model,v_l,T,x)
    p_v = pressure(model,v_v,T,y)
    μ_l = similar(F,n_c)
    μ_l = VT_chemical_potential!(μ_l,model,v_l,T,x)
    for i in 1:n_c
        F[i] = μ_l[i]
    end
    RT⁻¹ = 1/(Rgas(model)*Ts)
    μ_v = VT_chemical_potential!(μ_l,model,v_v,T,y)
    for i in 1:n_c
        μli = F[i]
        μvi = μ_v[i]
        Δμ = μli - μvi
        F[i] = Δμ*RT⁻¹
    end
    F[n_c+1] = (p_l-p_v)/ps
    return F
end

#non-condensable/non-volatile version
function μp_equality(model_long::EoSModel,model_short::EoSModel, F, T, v_long, v_short, x_long, x_short,ps_long,short_view, Ts = T)
    n_short = length(x_short)
    n_long = length(x_long)
    μ_long = similar(F,n_long)
    μ_long = VT_chemical_potential!(μ_long,model_long,v_long,T,x_long)
    p_long = pressure(model_long,v_long,T,x_long)
    p_short = pressure(model_short,v_short,T,x_short)
    RT⁻¹ = 1/(Rgas(model_long)*Ts)
    μ_long_view = @view(μ_long[short_view])
    for i in 1:n_short
        F[i] = μ_long_view[i]
    end
    μ_short = resize!(μ_long,n_short)
    μ_short = VT_chemical_potential!(μ_short,model_short,v_short,T,x_short)
    for i in 1:n_short
        μ_long_i = F[i]
        μ_short_i = μ_short[i]
        Δμ = (μ_long_i - μ_short_i)*RT⁻¹
        F[i] = Δμ
    end
    F[n_short+1] = (p_long-p_short)/ps_long
    return F
end

function μp_equality(model::EoSModel,::Nothing, F, T, v_l, v_v, x, y,ps,_view,Ts = T)
    return μp_equality(model,F,T,v_l,v_v,x,y,ps,Ts)
end

function VT_chemical_potential!(result,model,V,T,z)
    fun(x) = eos(model,V,T,x)
    return ForwardDiff.gradient!(result,fun,z)
end

function VT_chemical_potential_res!(result,model,V,T,z)
    fun(x) = eos_res(model,V,T,x)
    return ForwardDiff.gradient!(result,fun,z)
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

#include("general_eq")
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
