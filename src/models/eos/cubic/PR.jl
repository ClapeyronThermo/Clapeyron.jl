struct PRParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    acentricfactor::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type PRModel <: ABCubicModel end
@newmodel PR PRModel PRParam

export PR
function PR(components::Array{String,1}; userlocations::Array{String,1}=String[], idealmodel = BasicIdeal,verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    k  = params["k"]
    _pc = params["pc"]
    pc = _pc.values
    Tc = params["Tc"]
    Tc_ = Tc.values
    acentricfactor = params["w"]
    Ωa,Ωb = ab_consts(PRModel)
    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. Ωa*R̄^2*Tc_^2/pc), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. Ωb*R̄*Tc_/pc))

    packagedparams = PRParam(a, b, acentricfactor, Tc,_pc,Mw)
    model = PR(packagedparams,idealmodel)
    return model
end

function ab_consts(::Type{<:PRModel})
    return 0.457235,0.077796
end

function cubic_ab(model::PRModel,T,z=SA[1.0],n=sum(z))
    a = model.params.a.values
    b = model.params.b.values
    ω = model.params.acentricfactor.values
    Tc = model.params.Tc.values
    invn = one(n)/n
    αx = @. min((1+(0.37464+1.54226*ω-0.26992*ω^2)*(1-√(T/Tc)))^2,one(T)) * z * invn 
    āᾱ = dot(αx,Symmetric(a),αx)
    b̄ = dot(z,Symmetric(b),z) * invn*invn
    return āᾱ ,b̄
end

function cubic_abp(model::PRModel, V, T, z) 
    n = sum(z)
    v = V/n
    āᾱ ,b̄ = cubic_ab(model,T,z,n)
    _1 = one(b̄)
    denom = evalpoly(v,(-b̄*b̄,2*b̄,_1))
    p = R̄*T/(v-b̄) - āᾱ /denom
    return āᾱ, b̄,p
end

function cubic_poly(model::PRModel,p,T,z)
    a,b = cubic_ab(model,T,z)
    RT⁻¹ = 1/(R̄*T)
    A = a*p*RT⁻¹*RT⁻¹
    B = b*p*RT⁻¹
    k₀ = B*(B*(B+1.0)-A)
    k₁ = -B*(3*B+2.0) + A
    k₂ = B-1.0
    k₃ = 1.0
    return [k₀,k₁,k₂,k₃]
end

#=
 (-B2-2(B2+B)+A)
 (-B2-2B2-2B+A)
 (-3B2-2B+A)
=#
function a_res(model::PRModel, V, T, z)
    n = sum(z)
    a,b = cubic_ab(model,T,z,n)
    Δ1 = 1+√2
    Δ2 = 1-√2
    ΔPRΔ = 2*√2
    RT⁻¹ = 1/(R̄*T)
    ρ = n/V
    return -log(1-b*ρ) - a*RT⁻¹*log((Δ1*b*ρ+1)/(Δ2*b*ρ+1))/(ΔPRΔ*b)

    #return -log(V-n*b̄) + āᾱ/(R̄*T*b̄*2^(3/2)) * log((2*V-2^(3/2)*b̄*n+2*b̄*n)/(2*V+2^(3/2)*b̄*n+2*b̄*n))
end

cubic_zc(::PRModel) = 0.3074
