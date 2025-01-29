abstract type AlphaModel <:EoSModel end
abstract type TranslationModel <:EoSModel end

struct ABCubicParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct ABCCubicParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

const ONLY_VC = vcat(IGNORE_HEADERS,["Tc","Pc", "w"])
const ONLY_ACENTRICFACTOR = vcat(IGNORE_HEADERS,["Tc", "Pc", "Vc"])
"""
    ab_premixing(model,mixing,kij = nothing,lij = nothing)

given a model::CubicModel, that has `a::PairParam`, `b::PairParam`, a mixing::MixingRule and `kij`,`lij` matrices, `ab_premixing` will perform an implace calculation
to obtain the values of `a` and `b`, containing values aᵢⱼ and bᵢⱼ. by default, it performs the van der Wals One-Fluid mixing rule. that is:
```
aᵢⱼ = sqrt(aᵢ*aⱼ)*(1-kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
```
"""
function ab_premixing end

function ab_premixing(model::CubicModel,mixing::MixingRule,k = nothing, l = nothing)
    Ωa, Ωb = ab_consts(model)
    _Tc = model.params.Tc
    _pc = model.params.Pc
    a = model.params.a
    b = model.params.b
    diagvalues(a) .= @. Ωa*R̄^2*_Tc^2/_pc
    diagvalues(b) .= @. Ωb*R̄*_Tc/_pc
    epsilon_LorentzBerthelot!(a,k)
    sigma_LorentzBerthelot!(b,l)
    return a,b
end

function ab_premixing(model::CubicModel,kij::K,lij::L) where K <: Union{Nothing,PairParameter,AbstractMatrix} where L <: Union{Nothing,PairParameter,AbstractMatrix}
    return ab_premixing(model,model.mixing,kij,lij)
end

#legacy reasons
function ab_premixing(model::CubicModel,mixing::MixingRule,Tc,Pc,kij,lij)
    Ωa, Ωb = ab_consts(model)
    comps = Tc.components
    n = length(Tc)
    a = PairParam("a",comps,zeros(n,n))
    b = PairParam("b",comps,zeros(n,n))
    diagvalues(a) .= @. Ωa*R̄^2*Tc^2/pc
    diagvalues(b) .= @. Ωb*R̄*Tc/pc
    epsilon_LorentzBerthelot!(a,kij)
    sigma_LorentzBerthelot!(b,lij)
    return a,b
end

ab_premixing(model::CubicModel,mixing::MixingRule,Tc,pc,vc,kij,lij) = ab_premixing(model,mixing,Tc,pc,kij,lij) #ignores the Vc unless dispatch

function recombine_cubic!(model::CubicModel,k = nothing,l = nothing)
    recombine_mixing!(model,model.mixing,k,l)
    recombine_translation!(model,model.translation)
    recombine_alpha!(model,model.alpha)
    return model
end

function recombine_impl!(model::CubicModel)
    recombine_cubic!(model)
end

function c_premixing end

function cubic_ab(model::ABCubicModel,V,T,z=SA[1.0],n=sum(z))
    a = model.params.a.values
    b = model.params.b.values
    T = T * float(one(T))
    α = @f(α_function, model.alpha)
    c = @f(translation, model.translation)
    if length(z) > 1
        ā, b̄, c̄ = @f(mixing_rule, model.mixing, α, a, b, c)
    else
        ā = a[1, 1] * α[1]
        b̄ = b[1, 1]
        c̄ = c[1]
    end
    return ā, b̄, c̄
end

function data(model::ABCubicModel, V, T, z)
    n = sum(z)
    ā, b̄, c̄ = cubic_ab(model, V, T, z, n)
    return n, ā, b̄, c̄
end

function cubic_get_k end
function cubic_get_l end

get_k(model::CubicModel) = cubic_get_k(model,model.mixing,model.params)
get_l(model::CubicModel) = cubic_get_l(model,model.mixing,model.params)

function set_k!(model::CubicModel,k)
    check_arraysize(model,k)
    recombine_mixing!(model,model.mixing,k,nothing)
    return nothing
end

function set_l!(model::CubicModel,l)
    check_arraysize(model,k)
    recombine_mixing!(model,model.mixing,nothing,l)
    return nothing
end

function a_res(model::ABCubicModel, V, T, z,_data = data(model,V,T,z))
    n,ā,b̄,c̄ = _data
    Δ1,Δ2 = cubic_Δ(model,z)
    ΔΔ = Δ2 - Δ1
    RT⁻¹ = 1/(R̄*T)
    ρt = (V/n+c̄)^(-1) # translated density
    ρ  = n/V
    b̄ρt = b̄*ρt
    a₁ = -log1p((c̄-b̄)*ρ)
    if Δ1 == Δ2
        return a₁ - ā*ρt*RT⁻¹/(1-Δ1*b̄ρt)
    else
        l1 = log1p(-Δ1*b̄ρt)
        l2 = log1p(-Δ2*b̄ρt)
        return a₁ - ā*RT⁻¹*(l1-l2)/(ΔΔ*b̄)
    end
end

function cubic_poly(model::ABCubicModel,p,T,z)
    a,b,c = cubic_ab(model,p,T,z)
    RT⁻¹ = 1/(R̄*T)
    A = a*p*RT⁻¹*RT⁻¹
    B = b*p*RT⁻¹
    Δ1,Δ2 = cubic_Δ(model,z)
    ∑Δ = -Δ1 - Δ2
    Δ1Δ2 = Δ1*Δ2
    k₀ = -B*evalpoly(B,(A,Δ1Δ2,Δ1Δ2))
    k₁ = evalpoly(B,(A,-∑Δ,Δ1Δ2-∑Δ))
    k₂ = (∑Δ - 1)*B - 1
    k₃ = one(A) # important to enable autodiff
    return (k₀,k₁,k₂,k₃),c
end


function cubic_p(model::ABCubicModel, V, T, z,_data = @f(data))
    Δ1,Δ2 = cubic_Δ(model,z)
    n,a,b,c = _data
    v = V/n+c
    p = R̄*T/(v-b) - a/((v-Δ1*b)*(v-Δ2*b))
    return p
end

function pure_cubic_zc(model::ABCubicModel)
    Δ1,Δ2 = cubic_Δ(model,SA[1.0])
    _,Ωb = ab_consts(model)
    Ωb = only(Ωb)
    return (1 + (Δ1+Δ2+1)*Ωb)/3
end

function pure_cubic_zc(model::ABCCubicModel)
    Vc = model.params.Vc.values[1]
    pc = model.params.Pc.values[1]
    Tc = model.params.Tc.values[1]
    return pc*Vc/(R̄*Tc)
end

function second_virial_coefficient_impl(model::ABCubicModel,T,z = SA[1.0])
    a,b,c = cubic_ab(model,1/sqrt(eps(float(T))),T,z)
    return b - c - a/(Rgas(model)*T)
end

function lb_volume(model::CubicModel, z)
    V = 1e-5
    T = 0.0
    n = sum(z)
    invn = one(n) / n
    b = model.params.b.values
    c = @f(translation, model.translation)
    b̄ = dot(z, Symmetric(b), z) * invn #b has m3/mol units, result should have m3 units
    c̄ = dot(z, c) #result here should also be in m3
    return b̄ - c̄
end
#dont use αa, just a, to avoid temperature dependence
function T_scale(model::CubicModel, z)
    n = sum(z)
    invn2 = one(n) / (n * n)
    _Tc = model.params.Tc.values
    Tc = dot(z, _Tc) * invn2
    return Tc
end

function p_scale(model::CubicModel, z)
    n = sum(z)
    invn2 = one(n) / (n * n)
    _pc = model.params.Pc.values
    pc = dot(z, _pc) * invn2
    return pc
end

function x0_crit_pure(model::CubicModel)
    lb_v = lb_volume(model)
    (1.0, log10(lb_v / 0.3))
end

#works with models with a fixed (Tc,Pc) coordinate
function crit_pure_tp(model)
    single_component_check(crit_pure,model)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    Vc = volume(model,Pc,Tc,SA[1.])
    return (Tc,Pc,Vc)
end

function crit_pure_tp(model::ABCCubicModel)
    single_component_check(crit_pure,model)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    Vc = model.params.Vc.values[1]
    return (Tc,Pc,Vc)
end

function volume_impl(model::ABCubicModel,p,T,z,phase,threaded,vol0)
    lb_v = lb_volume(model,T,z)
    if iszero(p) && is_liquid(phase) #liquid root at zero pressure if available
        vl,_ = zero_pressure_impl(model,T,z)
        return vl
    elseif iszero(p) && is_vapour(phase) #zero pressure, ideal gas limit.
        _0 = zero(Base.promote_eltype(model,p,T,z))
        _1 = one(_0)
        return _1/_0
    end
    nRTp = sum(z)*R̄*T/p
    _poly,c̄ = cubic_poly(model,p,T,z)
    c = c̄*sum(z)
    num_isreal, z1, z2, z3 = Solvers.real_roots3(_poly)
    if num_isreal == 2
        vvl,vvg = nRTp*z1,nRTp*z2
    elseif num_isreal == 3
        vvl,vvg = nRTp*z1,nRTp*z3
    else
        vvl,vvg = nRTp*z1,nRTp*z1
    end
    #err() = @error("model $model Failed to converge to a volume root at pressure p = $p [Pa], T = $T [K] and compositions = $z")
    if !isfinite(vvl) && !isfinite(vvg) && phase != :unknown
        V0 = x0_volume(model, p, T, z; phase)
        v = _volume_compress(model, p, T, z, V0)
        #isnan(v) && err()
        return v
    end
    if num_isreal == 3 # 3 real roots
        vg = vvg
        _vl = vvl
        vl = ifelse(_vl > lb_v, _vl, vg) #catch case where solution is unphysical
    else # 1 real root (or 2 with the second one degenerate)
        vg = vl = z1 * nRTp
    end

    function gibbs(v)
        _df, f = ∂f(model, v, T, z)
        dv, dt = _df
        if abs((p + dv) / p) > 0.03
            return one(dv) / zero(dv)
        else
            return f + p * v
        end
    end
    #this catches the supercritical phase as well
    if vl ≈ vg
        return vl - c
    end

    if is_liquid(phase)
        return vl - c
    elseif is_vapour(phase)
        return vg - c
    else
        gg = gibbs(vg - c)
        gl = gibbs(vl - c)
        return ifelse(gg < gl, vg - c, vl - c)
    end
end

function pure_spinodal(model::ABCubicModel,T::K,v_lb::K,v_ub::K,phase::Symbol,retry,z = SA[1.0]) where K
    #=
    Segura, H., & Wisniak, J. (1997). Calculation of pure saturation properties using cubic equations of state. Computers & Chemical Engineering, 21(12), 1339–1347. doi:10.1016/s0098-1354(97)00016-1 
    =#
    a,b,c = cubic_ab(model,v_lb,T,z)
    Δ1,Δ2 = cubic_Δ(model,z)
    c1_c2 = - Δ1 - Δ2
    c1c2 = Δ1*Δ2
    RT = Rgas(model)*T
    bRT = b*RT
    Q4 = -RT
    Q3 = 2*(a - bRT*c1_c2)
    Q2 = b*(a*(c1_c2 - 4) - bRT*(c1_c2*c1_c2 + 2*c1c2))
    Q1 = 2*b*b*(a*(1 - c1_c2) - bRT*c1c2*c1_c2)
    Q0 = b*b*b*(a*c1_c2 - bRT*c1c2*c1c2)
    dpoly = (Q0,Q1,Q2,Q3,Q4)

    #on single component, a good approximate for vm is the critical volume.
    d2poly = (Q1,2*Q2,3*Q3,4*Q4)
    f = Base.Fix2(evalpoly,dpoly)
    nr,v1,v2,v3 = Solvers.real_roots3(d2poly)
    vroots = (v1,v2,v3)
    vm0 = findfirst(y -> (f(y) > 0 && y > b),vroots)
    if isnothing(vm0)
       return zero(v1)/zero(v1) 
    end
    vm = vroots[vm0]
    B = b - a/RT
    vx = ifelse(is_liquid(phase),b,-10B)
    v_bracket = minmax(vx,vm)
    prob = Roots.ZeroProblem(Base.Fix2(evalpoly,dpoly),v_bracket)
    vs = Roots.solve(prob)
    return vs - c
end

function liquid_spinodal_zero_limit(model::ABCubicModel,z)
    R̄ = Rgas(model)
    function F(Tx)
        a,b,c = cubic_ab(model,0,Tx,z)
        Δ1,Δ2 = cubic_Δ(model,z)
        Ax = R̄*Tx
        Bx = -(Ax*b*(Δ1+Δ2) + a)
        Cx = b*(Ax*Δ1*Δ2*b + a)
        return Bx^2 - 4*Ax*Cx
    end
    T0 = T_scale(model,z)
    prob = Roots.ZeroProblem(F,T0)
    T = Roots.solve(prob)
    _,vl = zero_pressure_impl(model,T,z)
    return T,vl
end

function zero_pressure_impl(T,a,b,c,Δ1,Δ2,z)
    #0 = R̄*T/(v-b) - a/((v-Δ1*b)*(v-Δ2*b))
    #f(v) = ((v-Δ1*b)*(v-Δ2*b))*R̄*T - (v-b)*a
    #RT(v^2 -(Δ1+Δ2)vb + Δ1Δ2b2) - av + ab
    #RTv^2 -(RT*Δ1b+Δ2b - a)*v + (RT*Δ1Δ2b2 + ab)
    A = R̄*T
    B = -(R̄*T*b*(Δ1+Δ2) + a)
    C = b*(R̄*T*Δ1*Δ2*b + a)
    #Δ = B2 - 4AC
    #R̄*T*b*(Δ1+Δ2)^2 + 2*R̄*T*b*(Δ1+Δ2)*a + a2 - 4*R̄*T*b*(R̄*T*Δ1*Δ2*b + a)
    #R̄*T*b*(Δ1+Δ2)^2 + 2*R̄*T*b*(Δ1+Δ2)*a + a2 - 4*R̄*T*b*(R̄*T*Δ1*Δ2*b + a)
    Δ = sqrt(B^2 - 4*A*C)
    vl = (-B - Δ)/(2*A) - c
    vmax = -B/(2*A) - c
    return vl,vmax
end

function zero_pressure_impl(model::ABCubicModel,T,z)
    a,b,c = cubic_ab(model,0,T,z)
    Δ1,Δ2 = cubic_Δ(model,z)
    return zero_pressure_impl(T,a,b,c,Δ1,Δ2,z)
end

function ab_consts(model::CubicModel)
    return ab_consts(typeof(model))
end

has_fast_crit_pure(model::ABCubicModel) = true

function x0_saturation_temperature(model::ABCubicModel,p,::Nothing)
    crit = crit_pure(model)
    return x0_saturation_temperature_crit(model, p, crit)
end

#=
#on the dpdv limit:
dp/dv = 0
p = RT/v-b - a/pol(v)
dpdv = -RT/(v-b)^2 + a/pol^2 * dpol = a*k -RT/(v-b)^2

vdw: pol = v2 -> pol(b) = b2, dpol(b) = 2b
pr: pol = v2 + 2bv - b2 -> pol(b) = 2b2, dpol(b) = 2v + 2b = 4b
rk: pol = v*(v+b) -> pol(b) = 2b2, dpol(b) = 2v + b = 3b

vdw:k = 2b/(b2)^2 = 2/b3 , k^-1 = 0.5b3
pr:k = 4b/(2b^2) = 1/b3, k^-1 = b3
rk:k = 3b/(2b^2) = 0.75/b3 lower  1.33b3

we want the lowest possible volume, to be sure on being on the liquid side.

solving for dpdv = 0
0 = a*k -RT/(v-b)^2
(v-b)^2 = RT/ak
v2 - 2vb + b2 - RT/ak = 0
v = b ± sqrt(b2 +  RT/ak - b2) #v > b
v = b + sqrt(kb3RT/a)
the lowest volume is reached with k(vdw):
vl = b + sqrt(0.5RTb3/2a)
on models with translation:
vl = b + sqrt(0.5RTb3/2a) - c
=#


function wilson_k_values!(K,model::ABCubicModel, p, T, crit = nothing)
    Pc = model.params.Pc.values
    Tc = model.params.Tc.values
    α = typeof(model.alpha)
    w1 = getparam(model,:acentricfactor)
    w2 = getparam(model.alpha,:acentricfactor)

    #we can find stored acentric factor values, so we calculate those
    if w1 !== nothing
        ω = w1.values
    elseif w2 !== nothing
        ω = w2.values
    else
        pure = split_model(model)
        ω = zero(Tc)
        for i in 1:length(Tc)
            ps = first(saturation_pressure(pure[i], 0.7 * Tc[i]))
            ω[i] = -log10(ps / Pc[i]) - 1.0
        end
    end

    return @.K .= Pc / p * exp(5.373 * (1 + ω) * (1 - Tc / T))

end

function vdw_tv_mix(Tc,Vc,z)
    Tm = zero(first(Tc)+first(Vc))
    Vm = zero(eltype(Vc))
    n = sum(z)
    invn2 = (1/n)^2
    for i in 1:length(z)
        zi = z[i]
        Vi = Vc[i]
        Ti = Tc[i]
        zii = zi*zi
        Vm += zii*Vi
        Tm += zii*Ti*Vi
        for j in 1:i-1
            zj = z[j]
            Vj = Vc[j]
            Tj = Tc[j]
            Tij = sqrt(Ti*Tj)
            Vij = 0.5*(Vi+Vj)
            zij = zj*zi
            Vm += 2zij*Vij
            Tm += zij*Tij*Viij
        end
    end
    Tcm = Tm/Vm
    Vcm = Vm*invn2
    return (Tcm,Vcm)
end

antoine_coef(model::ABCubicModel) = (6.668322465137264,6.098791871032391,-0.08318016317721941)
