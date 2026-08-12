abstract type AlphaModel <:EoSModel end
abstract type TranslationModel <:EoSModel end

Base.eltype(model::CubicModel) = Base.promote_eltype(model.alpha,model.translation,model.mixing,model.params)

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

struct SimpleAlphaParam <: EoSParam
    acentricfactor::SingleParam{Float64}
end

const ONLY_VC = vcat(IGNORE_HEADERS,["Tc","Pc", "w"])
const ONLY_ACENTRICFACTOR = vcat(IGNORE_HEADERS,["Tc", "Pc", "Vc"])
"""
    ab_premixing(model,mixing,kij = nothing,lij = nothing)

Given a model::CubicModel, that has `a::PairParam`, `b::PairParam`, a mixing::MixingRule and `kij`,`lij` matrices, `ab_premixing` will perform an implace calculation
to obtain the values of `a` and `b`, containing values aᵢⱼ and bᵢⱼ. By default, it performs the Van der Waals One-Fluid mixing rule. That is:
```
aᵢⱼ = sqrt(aᵢ*aⱼ)*(1-kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
```
"""
function ab_premixing end

"""
    ab_diagvalues!(a,b,Ωa,Ωb,Tc,Pc,R)
    ab_diagvalues!(model)

Calculates the diagonal (pure) terms of `a` and `b` in a cubic model, ignoring non-missing entries.
"""
function ab_diagvalues!(a::PairParam,b::PairParam,Ωa::Number,Ωb::Number,Tc,Pc,R̄)
    for i in 1:length(Tc)
        Tci,Pci = Tc[i],Pc[i]
        if a.ismissingvalues[i,i]
            a[i,i] = Ωa*R̄^2*Tci^2/Pci
        end

        if b.ismissingvalues[i,i]
            b[i,i] = Ωb*R̄*Tci/Pci
        end
    end
    return nothing
end

function ab_diagvalues!(a::PairParam,b::PairParam,Ωa::AbstractVector,Ωb::AbstractVector,Tc,Pc,R̄)
    for i in 1:length(Tc)
        Tci,Pci = Tc[i],Pc[i]
        if a.ismissingvalues[i,i]
            a[i,i] = Ωa[i]*R̄^2*Tci^2/Pci
        end

        if b.ismissingvalues[i,i]
            b[i,i] = Ωb[i]*R̄*Tci/Pci
        end
    end
    return nothing
end

function ab_diagvalues!(model::EoSModel)
    Ωa, Ωb = ab_consts(model)
    Tc = model.params.Tc
    Pc = model.params.Pc
    a = model.params.a
    b = model.params.b
    ab_diagvalues!(a,b,Ωa,Ωb,Tc,Pc,Rgas(model))
end

function ab_premixing(model::EoSModel,mixing::MixingRule, k, l)
    a = model.params.a
    b = model.params.b
    ab_diagvalues!(model)
    epsilon_LorentzBerthelot!(model.params.a,k)
    sigma_LorentzBerthelot!(model.params.b,l)
    return a,b
end

ab_premixing(model::EoSModel,mixing::MixingRule,k) = ab_premixing(model,mixing,k,nothing)
ab_premixing(model::EoSModel,mixing::MixingRule) = ab_premixing(model,mixing,nothing,nothing)

function ab_premixing(model::EoSModel,kij::K,lij::L) where K <: Union{Nothing,PairParameter,AbstractMatrix} where L <: Union{Nothing,PairParameter,AbstractMatrix}
    return ab_premixing(model,model.mixing,kij,lij)
end

function recombine_cubic!(model::CubicModel,k = nothing,l = nothing)
    recombine_mixing!(model,model.mixing,k,l)
    recombine_translation!(model,model.translation)
    recombine_alpha!(model,model.alpha)
    return model
end

function recombine_impl!(model::CubicModel)
    recombine_cubic!(model)
end

c_premixing(model) = nothing

function cubic_ab(model::CubicModel,V,T,z=SA[1.0])
    a = model.params.a.values
    b = model.params.b.values
    T = T * float(one(T))
    α = @f(α_function, model.alpha)

    if length(z) > 1
        return @f(mixing_rule, model.mixing, α, a, b)
    else
        return @f(mixing_rule1, model.mixing, α, a, b)
    end
end


function mixing_rule(model,V,T,z,mixing_model,α,a,b,c)
    c̄ = dot(z,c)/sum(z)
    ā,b̄,_ = @f(mixing_rule, model.mixing, α, a, b)
    return ā,b̄,c̄
end

#mixing rules: optimization for one-component
function mixing_rule1(model,V,T,z,mixing_model,α,a,b)
    _1 = oneunit(z[1])
    ā = a[1, 1] * α[1] * _1
    b̄ = b[1, 1] * _1
    c̄ = translation2(model,V/sum(z),T,SA[1.0],model.translation,ā,b̄,α)
    return ā, b̄, c̄
end

#For compatibility with earlier versions of Clapeyron, where we used to calculate translation as a vector
function translation2(model,V,T,z,translation_model,a,b,α)
    c = translation(model,V,T,z,translation_model)
    return dot(c,z)
end

function mixing_rule1(model,V,T,z,mixing_model,α,a,b,c)
    ā, b̄, _ = @f(mixing_rule1, model.mixing, α, a, b)
    c̄ = dot(c,z)/sum(z)
    return ā, b̄, c̄
end

function mixing_rule(model,V,T,z,mixing_model,α,a,b)
    c = @f(translation, model.translation)
    @f(mixing_rule, model.mixing, α, a, b, c)
end

function data(model::CubicModel, V, T, z)
    n = sum(z)
    ā, b̄, c̄ = cubic_ab(model, V, T, z)
    return Base.promote(n, ā, b̄, c̄)
end

get_k(model::CubicModel) = cubic_get_k(model,model.mixing,model.params)
get_l(model::CubicModel) = cubic_get_l(model,model.mixing,model.params)

cubic_get_k(model,mixing,params) = get_k_geomean(params.a.values)
cubic_get_l(model,mixing,params) = get_k_mean(params.b.values)

function set_k!(model::CubicModel,k)
    check_arraysize(model,k)
    recombine_mixing!(model,model.mixing,k,nothing)
    return nothing
end

function set_l!(model::CubicModel,l)
    check_arraysize(model,l)
    recombine_mixing!(model,model.mixing,nothing,l)
    return nothing
end

function a_res(model::DeltaCubicModel, V, T, z,_data = data(model,V,T,z))
    n,ā,b̄,c̄ = _data
    Δ1,Δ2 = cubic_ΔT(model,T,z)
    ΔΔ = Δ2 - Δ1
    RT⁻¹ = 1/(R̄*T)
    ρt = (V/n+c̄)^(-1) # translated density
    ρ  = n/V
    b̄ρt = b̄*ρt
    a₁ = -log1p((c̄-b̄)*ρ)
    if Δ1 == Δ2
        return a₁ - ā*ρt*RT⁻¹/(1-real(Δ1)*b̄ρt)
    else
        l1 = log1p(-Δ1*b̄ρt)
        l2 = log1p(-Δ2*b̄ρt)
        dl = l1 - l2
        return a₁ - ā*RT⁻¹*real(dl/(ΔΔ*b̄))
    end
end

function cubic_p(model::DeltaCubicModel, V, T, z,_data = @f(data),Δ = cubic_ΔT(model,T,z))
    Δ1,Δ2 = Δ
    n,a,b,c = _data
    v = V/n+c
    p = Rgas(model)*T/(v-b) - real(a/((v-Δ1*b)*(v-Δ2*b)))
    return p
end

function cubic_pure_zc(model::DeltaCubicModel)
    Tc = model.params.Tc[1]
    Pc = model.params.Pc[1]
    b = cubic_lb_volume(model,Tc,SA[1.0])
    Δ1,Δ2 = cubic_ΔT(model,Tc,SA[1.0])
    ∑Δ = real(Δ1 + Δ2)
    B = b*Pc/(Rgas(model)*Tc)
    return (1 + (∑Δ + 1)*B)/3 #Pc
end

function cubic_pure_zc(model::CubicModel)
    Tc = model.params.Tc[1]
    Pc = model.params.Pc[1]
    return volume(model,Pc,Tc,SA[1.0])
end

#=
function cubic_pure_zc(model::ABCubicModel)
    Δ1,Δ2 = cubic_Δ(model,SA[1.0])
    return cubic_pure_zc(Δ1,Δ2)
end

function cubic_pure_zc(Δ1::Number, Δ2::Number)
    r2m1 = 1.0 - Δ2
    r1m1 = 1.0 - Δ1
    t1 = cbrt(r1m1*r2m1*r2m1)
    t2 = cbrt(r2m1*r1m1*r1m1)
    ζc = (t1 + t2 + 1.0)
    x1 = (1.0 + Δ1 + Δ2)
    return ζc/(3.0*ζc - x1)
end =#

function second_virial_coefficient_impl(model::CubicModel,T,z = SA[1.0])
    a,b,c = cubic_ab(model,1/sqrt(eps(float(T))),T,z)
    return sum(z)*(b - c - a/(Rgas(model)*T))
end

function lb_volume(model::CubicModel,T,z)
    b̄ = cubic_lb_volume(model,T,z,model.mixing)
    c̄ = translation2(model,b̄,T,z,model.translation,nothing,b̄,nothing) #result here should also be in m3
    return b̄ - c̄
end

#some cubic mixing rules allow for T-dependent b.
#the default case is assume@s T-independency.
#the translation is added at the level of lb_volume
cubic_lb_volume(model,T,z) = cubic_lb_volume(model, T, z, model.mixing)

function cubic_lb_volume(model, T, z, mixing)
    V = 1e-5
    n = sum(z)
    invn = one(n) / n
    b = model.params.b.values
    b̄ = dot(z, Symmetric(b), z) * invn #b has m3/mol units, result should have m3 units
end
#dont use αa, just a, to avoid temperature dependence
function T_scale(model::CubicModel, z)
    _Tc = model.params.Tc.values
    return dot(z, _Tc) / sum(z)
end

function p_scale(model::CubicModel, z)
    _pc = model.params.Pc.values
    return dot(z, _pc) / sum(z)
end

function x0_crit_pure(model::CubicModel,z)
    Tc = T_scale(model,z)
    lb_v = lb_volume(model,Tc,z)/sum(z)
    (1.0, log10(lb_v / 0.3))
end

#by default, we assume Tc/Pc are fixed, Vc is variable.
function crit_pure(model::CubicModel)
    single_component_check(crit_pure,model)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    Vc = volume(model,Pc,Tc,SA[1.])
    return (Tc,Pc,Vc)
end

function crit_pure(model::DeltaCubicModel)
    single_component_check(crit_pure,model)

    if !has_fast_crit_pure(model)
        x0c = x0_crit_pure(model,SA[1.0])
        return crit_pure(model,x0c)
    end

    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    b = cubic_lb_volume(model,Tc,SA[1.0])
    Δ1,Δ2 = cubic_ΔT(model,Tc,SA[1.0])
    RT = Rgas(model)*Tc
    RTp = RT/Pc
    Vc0 = (RTp + (real(Δ1 + Δ2) + 1)*b)/3
    c = translation(model,Vc0,Tc,SA[1.0])
    Vc = Vc0 - c[1]
    #we know that in AB-cubics, the critical point is already determined.
    #model isa ABCubicModel && return (Tc,Pc,Vc)

    #for a general cubic model, we check if the critical pressure corresponds to the calculated pressure
    a = model.params.a[1,1]
    Pc_calculated = RT/(Vc0-b) - real(a/((Vc0-Δ1*b)*(Vc0-Δ2*b)))
    Pc_calculated ≈ Pc && return (Tc,Pc,Vc)

    #we failed. that means Pc is not the actual critical pressure. iterate (around Tc) and found Vc
    (Tc1,Pc1,Vc1) = __crit_pure_Δ(Tc,Vc,Rgas(model),a,b,Δ1,Δ2)
    if isnan(Pc1)
        return (Tc,Pc,Vc) #bail out
    end
    return (Tc1,Pc1,Vc1 - c[1])
end

#given fixed Tc, calculate Vc.
function __crit_pure_Δ(T,v0,R,a,b,Δ1,Δ2)
    f(_v) = __crit_pure_Δ_obj(T,_v,R,a,b,Δ1,Δ2)
    prob = Roots.ZeroProblem(f,v0)
    v = Roots.solve(prob,Roots.Newton())
    poly = real((v - Δ1*b)*(v - Δ2*b))
    p = R*T/(v - b) - a/poly
    return (T,p,v)
end

function __crit_pure_Δ_obj(T,v,R,a,b,Δ1,Δ2)
    RT = R*T
    poly = real((v - Δ1*b)*(v - Δ2*b))
    bb = real(-b*(Δ1 + Δ2))
    aRT = a/RT
    dpdv_scale = v*v/RT
    d2pdv2_scale = dpdv_scale*v
    dpoly = real((-b*(Δ1 + Δ2) + 2*v))
    dpdv = -RT/(v - b)^2 + a*dpoly/poly/poly
    d2pdv2 = 2RT/(v - b)^3 - 2a*(dpoly*dpoly/poly - 1)/(poly*poly)
    f = dpdv*dpdv_scale
    return dpdv*dpdv_scale,dpdv/d2pdv2
end

function volume_impl(model::DeltaCubicModel,p,T,z,phase,threaded,vol0)
    n,a,b,c = data(model,p,T,z)
    u,w = cubic_uwT(model,T,z)
    st,v1,v2 = cubic_poly_solver(a,b,p,Rgas(model),T,u,w,phase)
    C = n*c
    if st > 0 || st == -1
        vx = st == 1 ? v1 : v2
        return n*(vx - c)
    end
    v1 ≈ v2 && return n*(v1 - c)

    data2 = (n,a,b,zero(c))
    RT = Rgas(model)*T
    a1 = a_res(model,n*v1,T,z,data2)
    a2 = a_res(model,n*v2,T,z,data2)
    #if we ever add volume-dependent translations, they should be implemented here as -log1p(c̄2/V2)+log1p(c̄1/V1) instead of log(v2/v1)
    Δg = a1 - a2 + p*(v1 - v2)/RT + log(v2/v1)
    return Δg > 0 ? n*(v2 - c) : n*(v1 - c)
end

function cubic_poly_solver(a,b,p,R,T,u,w,phase)
    RT = R*T
    _pr = p/RT
    A = a*_pr/RT
    B = b*_pr
    _1 = one(Base.promote_eltype(A,B,u))
    ∅ = oftype(_1,NaN)
    pr = _pr*_1

    iszero(p) && !is_liquid(phase) && return 2,oftype(_1,Inf),oftype(_1,Inf) #fast handling of infinite vapour root

    #polynomial to calculate roots in η = b/V variable formulation.
    pη3 = fma_evalpoly(b,(_1*a,_1*RT*w,_1*p*w))# a + RT*b*w + p*b^2*w
    pη2 = -(a - RT*b*u + p*b^2*(w - u)) #fma_evalpoly(b,(_1*a,))
    pη1 = RT*b + p*b^2*(1-u)
    pη0 = -p*b^2
    poly_η = (pη0,pη1,pη2,pη3)
    ηc = -pη2/(3*pη3) #critical local η, we can decide if the root is liquid or gas.
    #WARNING: (this criteria fails with the anomalous second maxwell loop at high T)

    nr,η1,ηI,ηR,Δ = Solvers.__roots3(poly_η)
    
    #special case, 3 roots, but both the upper and lower roots are invalid. happens on Clapeyron.volume(EPPR78(["carbon dioxide"]), 3311.0e5, 400.0)
    if nr == 3 && η1 < 0 && ηR > 1 && (0 <= ηI <= ηR)
        nr,η1,ηI,ηR,Δ = 1,ηI,zero(ηI),zero(ηI),Δ
    end

    η_norm = maximum(abs,poly_η)
    Δ₀ = abs(Δ)/(η_norm*η_norm*η_norm)
    nr == 0 && return -1,∅,∅
    good_solve = Δ₀ > 1e-12

    st1,_vx1 = cubic_poly_solver_status(η1,ηc,phase),b/η1 #check status of root 1 (and calculate root 1)
    vx1 = st1 > 0 ? _vx1 : ∅
    good_solve && nr == 1 && (return max(0,st1),_vx1,_vx1) #no other root,both liquid an vapour roots converge to the same phase
    good_solve && nr == 2 && (return st1,vx1,vx1) #2 roots, return the single root unless the less stable root is requested

    st2 = cubic_poly_solver_status(ηR,ηc,phase) #check status of root 2. if there are two valid roots, then we asked for it or the gibbs criteria is needed
    vx2 = st2 > 0 ? b/ηR : ∅

    good_solve && st1 == -1 && (return st2,vx2,vx2) #if root 2 is requested, return root 2
    good_solve && st2 == -1 && (return st1,vx1,vx1) #if root 1 is requested, return root 1
    ηlo,ηhi = minmax(η1,ηR)
    if good_solve
        vl,vv = b/ηhi,b/ηlo
        return 0,vl,vv #use gibbs criterion to choose root
    end

    #polynomial to refine liquid root in volume.
    pv3 = p
    pv2 = -(p*(1 - u)*b + RT)
    pv1 = fma_evalpoly(b,(_1*a,- RT*_1*u,_1*p*(w - u)))
    pv0 = -b*fma_evalpoly(b,(a*_1,RT*w*_1,_1*p*w))
    poly_v = (pv0,pv1,pv2,pv3)

    #polynomial to refine liquid root in S = Z - 1
    ps3 = _1
    ps2 = 2 + (u - 1)*B
    ps1 = A + fma_evalpoly(B,(one(u),u - 2,w - u)) #1 + (u - 2)*B + (w - u)*B^2 + A
    ps0 = A*(1 - B) - B*fma_evalpoly(B,(one(u),u,w)) #A*(1 - B) - B*(1 + u*B + w*B^2)
    poly_s = (ps0,ps1,ps2,ps3)

    st11,_,vsol1 = cubic_poly_solver_refine(ηhi,ηc,poly_v,poly_s,pr,b,phase)
    v1 = st11 > 0 ? vsol1 : ∅
    nr == 1 && (return max(0,st11),vsol1,vsol1) #no other root,both liquid an vapour roots converge to the same phase
    nr == 2 && st1l > 0 && return (return st11,v1,v1) #2 roots, return the single root unless the less stable root is requested

    st22,_,vsol2 = cubic_poly_solver_refine(ηlo,ηc,poly_v,poly_s,pr,b,phase)
    v2 = st22 > 0 ? vsol2 : ∅
    st11 == -1 && (return st22,v2,v2) #if root 2 is requested, return root 2
    st22 == -1 && (return st11,v1,v1) #if root 1 is requested, return root 1
    vl,vv = minmax(v1,v2)
    return 0,vl,vv #use gibbs criterion to choose root
end

function cubic_poly_solver_status(η,ηc,phase,ignore_bounds = false)
    !isfinite(η) && return -1
    st0 = η > ηc ? 1 : 2
    !ignore_bounds && η > 1 && return -1
    !ignore_bounds && η < 0 && return -1
    is_liquid(phase) && st0 == 2 && return -1
    is_vapour(phase) && st0 == 1 && return -1
    return st0
end

function cubic_poly_solver_refine(η,ηc,poly_v,poly_s,pr,b,phase)
    st = cubic_poly_solver_status(η,ηc,phase,true)
    if st == -1
        nan = oftype(Y1,NaN)
        return st,nan,nan
    end
    if st == 1
        fl = Base.Fix2(fma_evalpoly,poly_v)
        prob_vl = Roots.ZeroProblem(fl,b/η)
        vx = Roots.solve(prob_vl)
        ηx = b/vx
    else
        fl = Base.Fix2(fma_evalpoly,poly_s)
        v0 = b/η
        S0 = v0*pr - 1
        prob_sv = Roots.ZeroProblem(fl,S0)
        SV = Roots.solve(prob_sv)
        vx = (SV + 1)/pr
        ηx = b/vx
    end
    new_st = cubic_poly_solver_status(ηx,ηc,phase)
    return new_st,ηx,vx
end

function pure_spinodal(model::DeltaCubicModel,T::K,v_lb::K,v_ub::K,phase::Symbol,retry,z = SA[1.0]) where K
    #=
    Segura, H., & Wisniak, J. (1997). Calculation of pure saturation properties using cubic equations of state. Computers & Chemical Engineering, 21(12), 1339–1347. doi:10.1016/s0098-1354(97)00016-1
    =#
    a,b,c = cubic_ab(model,v_lb,T,z)
    u,w = cubic_uwT(model,T,z)
    RT = Rgas(model)*T
    bRT = b*RT

    Q4 = -RT
    Q3 = 2*(a - bRT*u)
    Q2 = b*(a*(u - 4) - bRT*(u*u + 2*w))
    Q1 = 2*b*b*(a*(1 - u) - bRT*w*u)
    Q0 = b*b*b*(a*u - bRT*w*w)
    dpoly = (Q0,Q1,Q2,Q3,Q4)
    #on single component, a good approximate for vm is the critical volume.
    d2poly = (Q1,2*Q2,3*Q3,4*Q4)
    f = Base.Fix2(evalpoly,dpoly)
    nr,v1,v2,v3 = Solvers.real_roots3(d2poly)
    vroots = (v1,v2,v3)
    vm0 = nr == 1 ? nr : findfirst(y -> (f(y) > 0 && y > b),vroots)
    if isnothing(vm0)
       return zero(v1)/zero(v1)
    end
    vm = vroots[vm0]
    B = b - a/RT
    if vm < b
        return zero(v1)/zero(v1)
    end
    vx = ifelse(is_liquid(phase),b,-10B)
    v_bracket = minmax(vx,vm)
    prob = Roots.ZeroProblem(Base.Fix2(evalpoly,dpoly),v_bracket)
    vs = Roots.solve(prob)
    return sum(z)*(vs - c)
end

function liquid_spinodal_zero_limit(model::DeltaCubicModel,z)
    R̄ = Rgas(model)
    function F(Tx)
        a,b,c = cubic_ab(model,0,Tx,z)
        Δ1,Δ2 = cubic_ΔT(model,Tx,z)
        Ax = R̄*Tx
        Bx = -(Ax*b*(Δ1+Δ2) + a)
        Cx = b*(Ax*Δ1*Δ2*b + a)
        return real(Bx^2 - 4*Ax*Cx)
    end
    T0 = T_scale(model,z)
    prob = Roots.ZeroProblem(F,T0)
    T = Roots.solve(prob)
    _,vl = zero_pressure_impl(model,T,z)
    return T,vl
end

function zero_pressure_impl(model,T,z)
    return default_volume_impl(model,0.0,T,z,:liquid,false,nothing)
end

function zero_pressure_impl(model::DeltaCubicModel,T,z)
    a,b,c = cubic_ab(model,0,T,z)
    Δ1,Δ2 = cubic_ΔT(model,T,z)
    return zero_pressure_impl(T,a,b,c,Δ1,Δ2,z)
end

function zero_pressure_impl(T,a,b,c,Δ1,Δ2,z)
    #0 = R̄*T/(v-b) - a/((v-Δ1*b)*(v-Δ2*b))
    #f(v) = ((v-Δ1*b)*(v-Δ2*b))*R̄*T - (v-b)*a
    #RT(v^2 -(Δ1+Δ2)vb + Δ1Δ2b2) - av + ab
    #RTv^2 -(RT*Δ1b+Δ2b - a)*v + (RT*Δ1Δ2b2 + ab)
    RT = R̄*T
    A = one(RT)/b
    B = -((Δ1+Δ2) + a/(RT*b))
    C = (Δ1*Δ2*b + a/RT)
    #Δ = B2 - 4AC
    #R̄*T*b*(Δ1+Δ2)^2 + 2*R̄*T*b*(Δ1+Δ2)*a + a2 - 4*R̄*T*b*(R̄*T*Δ1*Δ2*b + a)
    #R̄*T*b*(Δ1+Δ2)^2 + 2*R̄*T*b*(Δ1+Δ2)*a + a2 - 4*R̄*T*b*(R̄*T*Δ1*Δ2*b + a)
    Δ = sqrt(B^2 - 4*A*C)
    vl = (-B - Δ)/(2*A) - c
    vmax = -B/(2*A) - c
    n = sum(z)
    return n*real(vl),n*real(vmax)
end

function x0_sat_pure_crit_info(model::ABCubicModel,T,crit,z = SA[1.0])
    #=
    saturation pressure approximation for cubics, near the critical point.

    implementation based on:
    Leibovici, C. F. (1993). Variant and invariant properties from cubic equations of state. Fluid Phase Equilibria, 84, 1–8. doi:10.1016/0378-3812(93)85114-2
    Sugie, H., Iwahori, Y., & Lu, B. C.-Y. (1989). On the application of cubic equations of state: Analytical expression for α/Tr and improved liquid density calculations. Fluid Phase Equilibria, 50(1–2), 1–20. doi:10.1016/0378-3812(89)80281-x
    =#
    RT = Rgas(model)*T
    Tc,Pc,_ = crit
    Δ1,Δ2 = cubic_ΔT(model,Tc,z)
    Ωa,Ωb = ab_consts(Δ1,Δ2)
    f = 1 - Ωa^(1/3)

    #=
    We use the Leibovici critical expansion, but we calculate Ψ as a function of Ωa, like sugie proposes.

    Sugie:
    α/Tr = 1 - f/(1 - f) * log(Pr/Tr)

    Leibovici:
    α/Tr = 1 + (1 - Ψ)/(2 + Ψ) * log(Tr/Pr)
    α/Tr = 1 + f/(1 - f) * log(Tr/Pr) #log(Tr/Pr) = -log(Pr/Tr)
    α/Tr = 1 + k * log(Tr/Pr)
    then:
    (1 - Ψ)/(2 + Ψ) = f/(1 - f)
    1 - Ψ - f + Ψf = 2f + Ψf
    Ψ = 1 - 3f
    =#

    k = real(f/(1 - f))
    Ψ = real(1 - 3*f)
    a0,_,_ = cubic_ab(model,Pc,Tc,z)
    at,b,c = cubic_ab(model,Pc,T,z)
    Tr = T/Tc
    α_over_Tr = at/(a0*Tr)
    lnTrPr = (α_over_Tr - 1)/k
    TrPr = exp(lnTrPr)
    PrTr = 1/TrPr
    Pr = Tr/TrPr
    p = Pc*Pr
    ε = sqrt(lnTrPr)
    Y0 = (1 - Ψ*PrTr)/3
    Y1 = sqrt((2 + Ψ)*(1 - Ψ))/3
    Y2 = evalpoly(Ψ,(4,1,1))/36
    Y3 = evalpoly(Ψ,(-48,64,87,6,-1))/(288*3*Y1)
    Y4 = evalpoly(Ψ,(-40,-44,-27,4,-1))/1296
    B = b*p/RT
    Yv = evalpoly(ε,(Y0,Y1,Y2,Y3,Y4))
    Yl = evalpoly(ε,(Y0,-Y1,Y2)) # evalpoly(ε,(Y0,-Y1,Y2,-Y3,Y4))
    Zl = Yl + B
    Zv = Yv + B
    n = sum(z)
    vv = n*(Zv*RT/p - c)
    vl = n*(Zl*RT/p - c)
    return p,vl,vv
end

#Δ1,Δ2 -> Ωa,Ωb infraestructure

#default: most models will use this

function cubic_ΔT(model,T,z)
    Δ1,Δ2 = cubic_Δ(model,z)
    return complex(Δ1),complex(Δ2)
end

function cubic_uwT(model,T,z)
    Δ1,Δ2 = cubic_ΔT(model,T,z)
    u = real(- Δ1 - Δ2)
    w = real(Δ1*Δ2)
    return u,w
end

#leibovici constants
function cubic_K(model,z)
    Δ1,Δ2 = cubic_Δ(model,z)
    u = - Δ1 - Δ2
    w = Δ1*Δ2
    return (1 + u + w)/(u + 2)^2
end

function cubic_Δ(model,z)
    return cubic_Δ(typeof(model))
end

cubic_Δ(model::EoSModel) = cubic_Δ(typeof(model))

function ab_consts(model::ABCubicModel,z)
    Δ1,Δ2 = cubic_Δ(model,z)
    return ab_consts(Δ1,Δ2)
end

function ab_consts(model::ABCubicModel)
    Δ1,Δ2 = cubic_Δ(model)
    return ab_consts(Δ1,Δ2)
end

Base.@assume_effects :foldable function ab_consts(::Type{T}) where T <: ABCubicModel
    Δ1,Δ2 = cubic_Δ(T)
    return ab_consts(Δ1,Δ2)
end

Base.@assume_effects :foldable function ab_consts(Δ1::Number, Δ2::Number)
    #calculate critical constants, from https://doi.org/10.1016/j.fluid.2012.05.008
    #code adapted from feos
    r2m1 = 1.0 - Δ2
    r1m1 = 1.0 - Δ1
    term1 = (r1m1*r2m1*r2m1)^(1/3)
    term2 = (r2m1*r1m1*r1m1)^(1/3)
    ζc = (term1 + term2 + 1.0)
    ηc = 1/ζc
    Ωb⁻¹ = 3.0*ζc - (1.0 + Δ1 + Δ2)
    Ωb2 = Ωb⁻¹*Ωb⁻¹
    Ωa = ζc*ζc*ζc*(1.0 - ηc*Δ1) * (1.0 - ηc*Δ2) * (2.0 - ηc*(Δ1 + Δ2)) /
        ((ζc - 1) * Ωb2)
    Ωb = 1/Ωb⁻¹
    return (Ωa, Ωb)
end

has_fast_crit_pure(model::DeltaCubicModel) = true

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

function wilson_k_values!(K,model::CubicModel, p, T, crit)
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
        pure = split_pure_model(model)
        ω = zero(Tc)
        for i in 1:length(Tc)
            ps = first(saturation_pressure(pure[i], 0.7 * Tc[i]))
            ω[i] = -log10(ps / Pc[i]) - 1.0
        end
    end

    return @.K .= Pc / p * exp(5.3726985503194395 * (1 + ω) * (1 - Tc / T))  #5.37 = log(10)*7/3

end

function tp_flash_fast_K0!(K,model::CubicModel,p,T,z)
    w1 = getparam(model,:acentricfactor)
    w2 = getparam(model.alpha,:acentricfactor)
    if w1 == nothing && w2 == nothing
        return false
    else
        wilson_k_values!(K,model, p, T, nothing)
        return true
    end
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

function x0_crit_mix(model::CubicModel,z)
    tci = model.params.Tc.values
    ∑z = sum(z)
    T_c  = prod(tci[i]^(z[i]/∑z) for i ∈ 1:length(model))
    P_c = dot(model.params.Pc.values,z)/∑z
    V_c = volume(model,P_c,T_c,z,phase = :v)/∑z
    return (log10(V_c),T_c)
end

#antoine_coef(model::ABCubicModel) = (6.668322465137264,6.098791871032391,-0.08318016317721941)

function transform_params(::Type{ABCubicParam},params,components)
    n = length(components)

    Tc = get!(params,"Tc") do
        SingleParam("Tc",components)
    end

    Pc = get!(params,"Pc") do
        SingleParam("Pc",components)
    end

    a = get!(params,"a") do
        aa = PairParam("a",components,zeros(Base.promote_eltype(Pc,Tc),n,n),ones(Bool,n,n))
    end
    a isa SingleParam && (params["a"] = PairParam(a))

    b = get!(params,"b") do
        PairParam("b",components,zeros(Base.promote_eltype(Pc,Tc),n,n),ones(Bool,n,n))
    end
    b isa SingleParam && (params["b"] = PairParam(b))

    Mw = get!(params,"Mw") do
        SingleParam("Mw",components)
    end
    return params
end

function transform_params(::Type{ABCCubicParam},params,components)
    n = length(components)
    transform_params(ABCubicParam,params,components)
    Tc = params["Tc"]
    Pc = params["Pc"]
    Vc = get!(params,"Vc") do
        SingleParam("Vc",components,zeros(Base.promote_eltype(Tc,Pc),n),fill(true,n))
    end

    c = get!(params,"c") do
        PairParam("c",components,zeros(Base.promote_eltype(Pc,Tc,Vc),n,n),ones(Bool,n,n))
    end
    c isa SingleParam && (params["c"] = PairParam(c))
    return params
end


"""
    CubicModel(cubicmodel::Type{T},params::Dict{String,ClapeyronParam},components;
    idealmodel = BasicIdeal,
    alpha = nothing,
    mixing = nothing,
    activity = nothing,
    translation = nothing,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false) where T <: CubicModel

## Input models
- `idealmodel`: Ideal Model
- `alpha`: Alpha model
- `mixing`: Mixing model
- `activity`: Activity Model, used in the creation of the mixing model.
- `translation`: Translation Model

## Description

Empty Cubic model constructor.
It requires specifying all model arguments.
"""
function CubicModel(cubicmodel::Type{T},params,components;
    idealmodel = BasicIdeal,
    alpha = nothing,
    mixing = nothing,
    activity = nothing,
    translation = nothing,
    userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false) where T <: CubicModel

    if CubicModel isa EoSModel
        return cubicmodel
    end

    _components = format_components(components)
    PARAM = parameterless_type(fieldtype(cubicmodel,:params))
    transform_params(PARAM,params,_components)
    transform_params(T,params,_components)
    init_mixing = init_model(mixing,components,activity,mixing_userlocations,activity_userlocations,verbose)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_alpha = init_alphamodel(alpha,components,params,alpha_userlocations,verbose)
    init_translation = init_model(translation,components,translation_userlocations,verbose)
    cubicparams = build_eosparam(PARAM,params)
    references = default_references(cubicmodel)
    return cubicmodel(_components,init_alpha,init_mixing,init_translation,cubicparams,init_idealmodel,references)
end

export CubicModel
