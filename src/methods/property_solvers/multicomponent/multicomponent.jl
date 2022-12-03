
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
    f(x) = eos(model,V,T,x)/(R̄*T)
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

function μp_equality(model::EoSModel, F, T, v_l, v_v, x, y,ts,ps)
    n_c = length(x)
    μ_l = similar(F,n_c)
    μ_l = VT_chemical_potential!(μ_l,model,v_l,T,x)
    for i in 1:n_c
        F[i] = μ_l[i]
    end
    μ_v = VT_chemical_potential!(μ_l,model,v_v,T,y)
    for i in 1:n_c
        μli = F[i]
        Δμ = (μli -μ_v[i])/(R̄*ts[i])
        F[i] = Δμ
    end
    p_l = pressure(model,v_l,T,x)
    p_v = pressure(model,v_v,T,y)
    F[n_c+1] = (p_l-p_v)/ps
    return F
end

#non-condensable/non-volatile version
function μp_equality(model_long::EoSModel,model_short::EoSModel, F, T, v_long, v_short, x_long, x_short,ts_short,ps_long,short_view)
    n_short = length(x_short)
    n_long = length(x_long)
    μ_long = similar(F,n_long)
    μ_long = VT_chemical_potential!(μ_long,model_long,v_long,T,x_long)
    μ_long_view = @view(μ_long[short_view])
    for i in 1:n_short
        F[i] = μ_long_view[i]
    end
    μ_short = resize!(μ_long,n_short)
    μ_short = VT_chemical_potential!(μ_short,model_short,v_short,T,x_short)
    for i in 1:n_short
        μlong_i = F[i]
        Δμ = (μlong_i -μ_short[i])/(R̄*ts_short[i])
        F[i] = Δμ
    end
    p_long = pressure(model_long,v_long,T,x_long)
    p_short = pressure(model_short,v_short,T,x_short)
    F[n_short+1] = (p_long-p_short)/ps_long
    return F
end

function μp_equality(model::EoSModel,::Nothing, F, T, v_l, v_v, x, y,ts,ps,_view)
    return μp_equality(model,F,T,v_l,v_v,x,y,ts,ps)
end

function VT_chemical_potential!(result,model,V,T,z)
    fun(x) = eos(model,V,T,x)
    return ForwardDiff.gradient!(result,fun,z)
end

function VT_chemical_potential_res!(result,model,V,T,z)
    fun(x) = eos_res(model,V,T,x)
    return ForwardDiff.gradient!(result,fun,z)
end

function PV_critical_temperature(model,p)
    Tc,_,vc = crit_pure(model)
    g(T) = p - pressure(model,vc,T)
    gi = Roots.ZeroProblem(g,Tc)
    T = Roots.solve(gi)
    return T
end

#returns saturation temperature if below crit_pure, if not, it returns
function _sat_Ti(model,p)
    pure = split_model(model)
    n = length(pure)
    Tsat = first.(saturation_temperature.(pure,p))
    for i ∈ 1:n
        if isnan(Tsat[i])
            T = PV_critical_temperature(pure[i],p)
            Tsat[i] = T
        end
    end
    return Tsat
end

function aprox_psat(pure,T,crit)
    coeff  = antoine_coef(pure)
    if coeff !== nothing
        A,B,C = coeff
    else
        A,B,C = (6.668322465137264,6.098791871032391,-0.08318016317721941)
    end
    Tc,Pc,Vc = crit
    T̄ = T/Tc
    if T > Tc
        return pressure(pure,Vc,Tc)
    else
        return exp(A-B/(T̄+C))*Pc
    end
end

function wilson_k_values(model::EoSModel,p,T,crit = nothing)
    n = length(model)
    pure = split_model.(model)
    if crit === nothing
        crit = crit_pure.(pure)
    end
    K0 = zeros(typeof(p+T),n)
    for i ∈ 1:n
        pure_i = pure[i]
        Tc,pc,_ = crit[i]
        ps = first(saturation_pressure(pure_i,0.7*Tc))
        ω = -log10(ps/pc) - 1.0
        K0[i] = exp(log(pc/p)+5.373*(1+ω)*(1-Tc/T))
    end
    return K0
end

function bubbledew_check(vl,vv,zin,zout)
    (isapprox(vl,vv) && isapprox(zin,zout)) && return false 
    !all(isfinite,zout) && return false
    !all(isfinite,(vl,vv)) && return false
    return true
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

export bubble_pressure_fug, bubble_temperature_fug, dew_temperature_fug, dew_pressure_fug

export bubble_pressure,    dew_pressure,    LLE_pressure,    azeotrope_pressure, VLLE_pressure
export bubble_temperature, dew_temperature, LLE_temperature, azeotrope_temperature, VLLE_temperature
export dew_pressure_fug_condensable, dew_temperature_fug_condensable
export bubble_pressure_fug_volatile, bubble_temperature_fug_volatile
export crit_mix, UCEP_mix, UCST_mix
