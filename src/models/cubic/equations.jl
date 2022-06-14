struct ABCubicParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

"""
    ab_premixing(::Type{T},mixing,Tc,pc,kij) where T <: ABCubicModel

given `Tc::SingleParam`, `pc::SingleParam`, `kij::PairParam` and `mixing <: MixingRule`, it will return 
`PairParam`s `a` and `b`, containing values aᵢⱼ and bᵢⱼ. by default, it performs the van der Wals One-Fluid mixing rule. that is:
```
aᵢⱼ = sqrt(aᵢ*aⱼ)*(1-kᵢⱼ)
bᵢⱼ = (bᵢ + bⱼ)/2
```
"""
function ab_premixing end

function ab_premixing(::Type{T},mixing,Tc,pc,kij) where T <: ABCubicModel
    Ωa, Ωb = ab_consts(T)
    _Tc = Tc.values
    _pc = pc.values
    a = epsilon_LorentzBerthelot(SingleParam(pc, @. Ωa*R̄^2*_Tc^2/_pc),kij)
    b = sigma_LorentzBerthelot(SingleParam(pc, @. Ωb*R̄*_Tc/_pc))
    return a,b
end

function cubic_ab(model::ABCubicModel,V,T,z=SA[1.0],n=sum(z))
    a = model.params.a.values
    b = model.params.b.values
    T = T*float(one(T))
    α = @f(α_function,model.alpha)
    c = @f(translation,model.translation)
    if length(z)>1
    ā,b̄,c̄ = @f(mixing_rule,model.mixing,α,a,b,c)
    else
         ā = a[1,1]*α[1]
         b̄ = b[1,1]
         c̄ = c[1]
    end
    return ā ,b̄, c̄
end

function data(model::ABCubicModel,V,T,z)
    n = sum(z)
    ā ,b̄, c̄ = cubic_ab(model,V,T,z,n)
    return n, ā ,b̄, c̄
end

function second_virial_coefficient(model::ABCubicModel,T::Real,z = SA[1.0])
    a,b,c = cubic_ab(model,1/sqrt(eps(float(T))),T,z)
    return b-a/(R̄*T)
end

function lb_volume(model::CubicModel,z = SA[1.0])
    V = 1e-5
    T = 0.
    n = sum(z)
    invn = one(n)/n
    b = model.params.b.values
    c = @f(translation,model.translation)
    b̄ = dot(z,Symmetric(b),z)*invn #b has m3/mol units, result should have m3 units
    c̄ = dot(z,c)*invn
    return b̄-c̄
end
#dont use αa, just a, to avoid temperature dependence
function T_scale(model::CubicModel,z=SA[1.0])
    n = sum(z)
    invn2 = one(n)/(n*n)
    Ωa,Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = dot(z, Symmetric(_a), z)*invn2/Ωa
    b = dot(z, Symmetric(_b), z)*invn2/Ωb
    return a/b/R̄
end

function p_scale(model::CubicModel,z=SA[1.0])
    n = sum(z)
    invn2 = (1/n)^2
    Ωa,Ωb = ab_consts(model)
    _a = model.params.a.values
    _b = model.params.b.values
    a = invn2*dot(z, Symmetric(_a), z)/Ωa
    b = invn2*dot(z, Symmetric(_b), z)/Ωb
    return a/ (b^2) # Pc mean
end

function x0_crit_pure(model::CubicModel)
    lb_v = lb_volume(model)
    (1.0, log10(lb_v/0.3))
end

function crit_pure(model::ABCubicModel)
    Tc = model.params.Tc.values[1]
    Pc = model.params.Pc.values[1]
    Vc = volume(model,Pc,Tc,SA[1.],phase=:v)
    return (Tc,Pc,Vc)
end

function volume_impl(model::ABCubicModel,p,T,z=SA[1.0],phase=:unknown,threaded=false)
    lb_v   =lb_volume(model,z)
    RTp = R̄*T/p
    _poly,c̄ = cubic_poly(model,p,T,z)
    sols = Solvers.roots3(_poly)
    function imagfilter(x)
        absx = abs(imag(x)) 
        return absx < 8*eps(typeof(absx))
    end
    x1,x2,x3 = sols
    sols = (x1,x2,x3)
    xx = (x1,x2,x3)
    isreal = imagfilter.(xx)
    vvv = extrema(real.(xx))
    
    zl,zg = vvv
    vvl,vvg = RTp*zl,RTp*zg
    err() = @error("model $model Failed to converge to a volume root at pressure p = $p [Pa], T = $T [K] and compositions = $z")
    if sum(isreal) == 3 #3 roots
        vg = vvg
        _vl = vvl
        vl = ifelse(_vl>lb_v,_vl,vg) #catch case where solution is unphysical
    elseif  sum(isreal) == 1
        i = findfirst(imagfilter,sols)
        vl = real(sols[i])*RTp
        vg = real(sols[i])*RTp
    elseif  sum(isreal) == 0
       
        V0 = x0_volume(model,p,T,z;phase)
        v = _volume_compress(model,p,T,z,V0)
        isnan(v) && err()
        return v
    end

    function gibbs(v)
        _df,f =  ∂f(model,v,T,z)
        dv,dt = _df
        if abs((p+dv)/p) > 0.03
            return one(dv)/zero(dv)
        else
            return f + p*v
        end
    end
    #this catches the supercritical phase as well
    if vl ≈ vg
        return vl-c̄
    end

    if is_liquid(phase)
        return vl-c̄
    elseif is_vapour(phase)
        return vg-c̄
    else
        gg = gibbs(vg-c̄)
        gl = gibbs(vl-c̄)
        return ifelse(gg<gl,vg-c̄,vl-c̄)
    end
end

function ab_consts(model::CubicModel)
    return ab_consts(typeof(model))
end

function x0_sat_pure(model::ABCubicModel,T)
    a,b,c = cubic_ab(model,1/sqrt(eps(float(T))),T)
    Tc = model.params.Tc.values[1]
    pc = model.params.Pc.values[1]
    zc = cubic_zc(model)
    vc = zc*R̄*Tc/pc - c
    if Tc < T
        nan = zero(T)/zero(T)
        return (nan,nan)
    end
    B = b-a/(R̄*T)
    pv0 = -0.25*R̄*T/B
    vl = b + sqrt(0.5R̄*T*b^3/a) - c
    pc = model.params.Pc.values[1]
    p_vl = pressure(model,vl,T)
    p_low = min(p_vl,pc)
    pl0 = max(zero(b),p_low)
    p0 = 0.5*(pl0+pv0)
    vv = volume_virial(B,p0,T) - c
    if p_vl > pc #improves predictions around critical point
        vlc,vvc =  vdw_x0_xat_pure(T,Tc,pc,vc)
        vl = 0.5*(vl+vlc)
        vv = 0.5*(vv+vvc)
    end
    return (log10(vl),log10(vv))
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
pr:k =  4b/(2b^2) = 1/b3, k^-1 = b3
rk:k =  3b/(2b^2) = 0.75/b3 lower  1.33b3

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

function wilson_k_values(model::ABCubicModel,p,T)
    Pc = model.params.Pc.values
    Tc = model.params.Tc.values

    if hasfield(typeof(model.alpha.params),:acentricfactor)
        ω = model.alpha.params.acentricfactor.values
    else
        pure = split_model(model)
        ω = zero(Tc)
        for i in 1:length(Tc)
            ps = first(saturation_pressure(pure[i],0.7*Tc[i]))
            ω[i] = -log10(ps/Pc[i]) - 1.0
        end
    end

    return  @. Pc/p*exp(5.373*(1+ω)*(1-Tc/T))
end
