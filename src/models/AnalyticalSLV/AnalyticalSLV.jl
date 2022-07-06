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

struct AnalyticalSLV{T <: IdealModel} <:AnalyticalSLVModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    params::AnalyticalSLVParam
    idealmodel::T
    references::Array{String,1}
end

@registermodel AnalyticalSLV

function AnalyticalSLV(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false)
    params = getparams(components, ["properties/critical.csv", 
                                    "properties/molarmass.csv",
                                    "AnalyticalSLV/SLV_Like.csv"];
                                    userlocations=userlocations,
                                    verbose=verbose)
    
    Mw = params["Mw"]
    Tc = params["Tc"]
    Pc = params["pc"]
    Vc = params["vc"]

    a0 = params["a0"]
    a1 = params["a1"]
    a2 = params["a2"]
    n =  params["n"]

    b0 = params["b0"]
    b1 = params["b1"]
    b2 = params["b2"]
    m =  params["m"]
    c = params["cr"]
    d = params["dr"]

    c.values .= c.values .* Vc.values
    d.values .= d.values .* Vc.values

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    packagedparams = AnalyticalSLVParam(Mw,Tc,Pc,Vc,a0,a1,a2,n,b0,b1,b2,m,c,d)
    references = ["10.1023/a:1024015729095"]
    icomponents = 1:length(components)
    model = AnalyticalSLV(components,icomponents,packagedparams,init_idealmodel,references)
    return model
end

function abcd(model::AnalyticalSLVModel,V,T,z=SA[1.0])
    n = sum(z)
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
    ā = ā/n/n
    b̄ = b̄/n
    c̄ = dot(c,z)/n
    d̄ = dot(d,z)/n
    return ā,b̄,c̄,d̄
end

function data(model::AnalyticalSLVModel,V,T,z)
    return abcd(model,V,T,z)
end

function a_res(model::AnalyticalSLVModel,V,T,z,_data = @f(data))
    ā,b̄,c̄,d̄ = _data
    n = sum(z)
    RT⁻¹ = 1/(R̄*T)
    ρ = n/V
    
    return ā*ρ*RT⁻¹
end

function x0_volume_liquid(model::AnalyticalSLVModel,T,z = SA[1.0])
    return 1.01*dot(model.params.c.values,z)
end

function x0_volume_solid(model::AnalyticalSLVModel,T,z = SA[1.0])
    b0 = model.params.b0.values
    b1 = model.params.b1.values
    b2 = model.params.b2.values
    m = model.params.m.values
    Tc = model.params.Tc.values
    Vc = model.params.Vc.values

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
    return 1.01*b̄
end

#(b*(b-d)*log(v-b) - c*(c-d)*log(V-c))/(b-c) - b + c - d 
