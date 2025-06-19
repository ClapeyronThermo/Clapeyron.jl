abstract type Kiselev2000Model <: CrossOverModel end

struct Kiselev2000Param{T} <: ParametricEoSParam{T}
    Tc::SingleParam{T}
    Vc::SingleParam{T}
    d1::SingleParam{T}
    v1::SingleParam{T}
    Gi::SingleParam{T}
end

function Z_base(model, V, T, z)
    ares(x)  = a_res(model, x, T, z)
    dares(x) = Solvers.derivative(ares,x)
    return 1.0 - V*dares(V)
end

@newmodelsimple Kiselev2000 Kiselev2000Model Kiselev2000Param

function a_res_crossover(model::CrossOver,V,T,z,critmodel::Kiselev2000Model)
    ∑z = sum(z)
    Tc = critmodel.params.Tc.values
    vc = critmodel.params.Vc.values
    Tc0 = model.params.Tc0.values
    vc0 = model.params.Vc0.values
    d1 = critmodel.params.d1.values
    v1 = critmodel.params.v1.values
    Gi = critmodel.params.Gi.values
    basemodel = model.basemodel
    _0 = zero(Base.promote_eltype(critmodel,basemodel,V,T,z))
    Δτ = _0
    Δϕ = _0
    _Tc0 = _0
    _vc0 = _0
    _d1 = _0
    _v1 = _0
    invGi = _0
    for i ∈ @comps
        zᵢ,Tcᵢ,vcᵢ,Tc0ᵢ,vc0ᵢ,d1ᵢ,v1ᵢ,Giᵢ = z[i],Tc[i],vc[i],Tc0[i],vc0[i],d1[i],v1[i],Gi[i]
        Δτ += zᵢ*(Tcᵢ/Tc0ᵢ - 1.0)
        Δϕ += zᵢ*(vcᵢ/vc0ᵢ - 1.0)
        _Tc0 += zᵢ*Tc0ᵢ
        _vc0 += zᵢ*vc0ᵢ
        _d1 += zᵢ*d1ᵢ
        _v1 += zᵢ*v1ᵢ
        invGi += zᵢ/Giᵢ
    end
    Δτ /= ∑z
    Δϕ /= ∑z
    _Tc0 /= ∑z
    _vc0 /= ∑z
    _Tc = _Tc0 + Δτ*_Tc0
    _vc = _vc0 + Δϕ*_vc0
    _d1 /= ∑z
    _v1 /= ∑z
    _Gi = ∑z/invGi
    b = sqrt(1.359)
    m = 1.3
    beta = 0.325
    alpha = 0.11
    gamma = 2.0-alpha-2.0*beta
    D = 0.51
    v = V/∑z
    τ = T/_Tc - 1.0
    ϕ = v/_vc - 1.0
    fnc = 4.0*(b/m*abs(ϕ*(1.0 + _v1*(ϕ^2)*exp(-8.5*ϕ)) + _d1*τ*(1.0 - 2.0*τ)))^(1.0/beta) + 2.0*τ
    q = sqrt((fnc + sqrt(fnc^2 + 12.0*τ^2))/6.0/_Gi)
    Y = (q/(1.0+q))^2
    τY = 0.0
    if Y > 0.0
        τY = τ*Y^(-alpha/2.0)
    end
    τs = τY + (1.0 + τ)*Δτ*Y^(2.0*(2.0 - alpha)/3.0)
    Ts = _Tc0*(τs + 1.0)
    ϕs = ϕ*Y^((gamma - 2.0*beta)/4.0) + (1.0 + ϕ)*Δϕ*Y^((2.0 - alpha)/2.0)
    vs = _vc0*(ϕs + 1.0)
    Δv = v/_vc0 - 1.0
    _a_res = a_res(basemodel, vs, Ts, z) - a_res(basemodel, _vc0, Ts, z) + a_res(basemodel, _vc0, T, z) + ϕs*Z_base(basemodel, _vc0, Ts, z) - Δv*Z_base(basemodel, _vc0, T, z) + log((Δv + 1.0)/(ϕs + 1.0))
    return _a_res
end

default_references(::Type{Kiselev2000}) = ["10.1023/A:1006657410862"]