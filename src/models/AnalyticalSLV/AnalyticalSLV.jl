struct AnalyticalSLVParam <: EoSParam
    Mw::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    a0::SingleParam{Float64}
    a1::SingleParam{Float64}
    a2::SingleParam{Float64}
    n::SingleParam{Float64}
    b0::SingleParam{Float64}
    b1::SingleParam{Float64}
    b2::SingleParam{Float64}
    m::SingleParam{Float64}
    c::SingleParam{Float64}
    d::SingleParam{Float64}
end

abstract type AnalyticalSLVModel <: EoSModel end
@newmodel AnalyticalSLV AnalyticalSLVModel AnalyticalSLVParam false
default_locations(::Type{AnalyticalSLV}) = ["properties/critical.csv","properties/molarmass.csv","AnalyticalSLV/SLV_Like.csv"]
default_references(::Type{AnalyticalSLV}) = ["10.1023/a:1024015729095"]
function transform_params(::Type{AnalyticalSLV},params)
    c = params["cr"]
    d = params["dr"]
    Vc = params["Vc"]
    c.values .= c.values .* Vc.values
    d.values .= d.values .* Vc.values
    params["c"] = c
    params["d"] = d
    return params
end

export AnalyticalSLV

function abcd(model::AnalyticalSLVModel,V,T,z=SA[1.0])
    ∑z = sum(z)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = model.params.Vc.values

    a0 = model.params.a0.values
    a1 = model.params.a1.values
    a2 = model.params.a2.values
    n = model.params.n.values

    b0 = model.params.b0.values
    b1 = model.params.b1.values
    b2 = model.params.b2.values
    m = model.params.m.values

    c = model.params.c.values
    d = model.params.d.values

    _0 = zero(T+first(z))

    ā = _0
    b̄ = _0
    for i in @comps
        Tci = Tc[i]
        Tri = T/Tci
        zi = z[i]
        br = b0[i] + b1[i]*exp(-b2[i]*Tri^m[i])
        bi = br*Vc[i]
        b̄ += bi*zi
        ari = a0[i] + a1[i]*exp(-a2[i]*Tri^n[i])
        ai = ari*(R̄*Tci)^2/Pc[i]
        ā += zi*zi*ai
        for j in 1:(i-1)
            Tcj = Tc[j]
            Trj = T/Tcj
            arj = a0[j] + a1[j]*exp(-a2[j]*Trj^n[j])
            aj = arj*(R̄*Tcj)^2/Pc[j]
            aij = sqrt(ai*aj)
            ā += zi*z[j]*aij
        end
    end
    ā = ā/∑z/∑z
    b̄ = b̄/∑z
    c̄ = dot(c,z)/∑z
    d̄ = dot(d,z)/∑z
    return ā,b̄,c̄,d̄
end

function x0_volume_solid(model::AnalyticalSLVModel,T,z = SA[1.0])
    ∑z = sum(z)
    Tc = model.params.Tc.values
    Vc = model.params.Vc.values

    b0 = model.params.b0.values
    b1 = model.params.b1.values
    b2 = model.params.b2.values
    m = model.params.m.values

    _0 = zero(T+first(z))

    b̄ = _0
    for i in @comps
        Tci = Tc[i]
        Tri = T/Tci
        zi = z[i]
        br = b0[i] + b1[i]*exp(-b2[i]*Tri^m[i])
        bi = br*Vc[i]
        b̄ += bi*zi
    end
    b̄ = b̄/∑z
    return 1.01*b̄
end

function lb_volume(model::AnalyticalSLVModel,z)
    return x0_volume_solid(model,0,z)/1.01
end

function data(model::AnalyticalSLVModel,V,T,z)
    return abcd(model,V,T,z)
end

function x0_volume_liquid(model::AnalyticalSLVModel,T,z = SA[1.0])
    return 1.01*dot(model.params.c.values,z)
end

function T_scale(model::AnalyticalSLVModel, z)
    n = sum(z)
    invn2 = one(n) / (n * n)
    _Tc = model.params.Tc.values
    Tc = dot(z, _Tc) * invn2
    return Tc
end

function x0_crit_pure(model::AnalyticalSLVModel)
    Pc = model.params.Pc.values[1]
    Tc = model.params.Tc.values[1]
    vc = volume(model,Pc,Tc,phase = :v)
    return (1.01,log10(2*vc))
end

function p_scale(model::AnalyticalSLVModel, z)
    n = sum(z)
    invn2 = one(n) / (n * n)
    _pc = model.params.Pc.values
    pc = dot(z, _pc) * invn2
    return pc
end

function _pressure(model::AnalyticalSLVModel,V,T,z=SA[1.0])
    _data = @f(data)
    ā,b̄,c̄,d̄ = _data
    v = V/sum(z)
    return R̄*T*(v-d̄)/(v-b̄)/(v-c̄) - ā/(v^2)
end

cubic_Δ(model::AnalyticalSLV,z) = (0.0,0.0)

function a_res(model::AnalyticalSLVModel,V,T,z,_data = @f(data))
    ā,b̄,c̄,d̄ = _data
    n = sum(z)
    RT⁻¹ = 1/(R̄*T)

    ρ = n/V
    Δ1,Δ2 = cubic_Δ(model,z)
    ΔΔ = Δ2 - Δ1

    bd = (b̄ - d̄)
    dc = (d̄ - c̄)
    k1 = bd*log(1-b̄*ρ)
    k2 = dc*log(abs(1-c̄*ρ))
    b̄ρt = b̄*ρ

    #The integral of 1/x is log(abs(x))
    #On solid volumes, 1-c̄*ρ is negative, so the abs matters
    a₁ = -(k1 + k2)/(b̄ - c̄) - (d̄ - c̄)/(b̄ - c̄)
    if Δ1 == Δ2
        return a₁ - ā*ρt*RT⁻¹/(1-Δ1*b̄ρt)
    else
        l1 = log1p(-Δ1*b̄ρt)
        l2 = log1p(-Δ2*b̄ρt)
        return a₁ - ā*RT⁻¹*(l1-l2)/(ΔΔ*b̄)
    end
   #return -(k1 + k2)/(b̄ - c̄) -  ā*ρ*RT⁻¹ - (d̄ - c̄)/(b̄ - c̄)
end

#a_res(V,T,z) = -integral(V -> (z(V)-1)/V,Inf,V)
#=

Px = RT/(V-b) * (V-d)/(V-c)
Zx = Tr(Vr - dr)/(Vr - br)/(Vr - cr)
integral(z/V) = (Tr*Vc*(br - dr)*(log(V - br*Vc)) + (dr -cr)*log(V -cr*Vc))/(br - cr)


=#