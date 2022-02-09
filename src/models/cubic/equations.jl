

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
    invn2 = (one(n)/n)^2
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
        c̄ = c[1,1]
    end
    return ā ,b̄, c̄
end

function second_virial_coefficient(model::ABCubicModel,T::Real,z = SA[1.0])
    #@info "fast shortcut"
    a,b,c = cubic_ab(model,1e-4,T,z)
    return b-a/(R̄*T)
end

function lb_volume(model::CubicModel,z = SA[1.0])
    V = 1e-5
    T = 0.
    n = sum(z)
    invn = one(n)/n
    b = model.params.b.values
    c = @f(translation,model.translation)
    b̄ = dot(z,Symmetric(b),z)*invn*invn
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

function volume(model::ABCubicModel,p,T,z=SA[1.0];phase=:unknown,threaded=false)
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


    


