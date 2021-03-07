struct PRParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    acentricfactor::SingleParam{Float64}
    Tc::SingleParam{Float64}
end

abstract type PRModel <: ABCubicModel end
@newmodel PR PRModel PRParam

export PR
function PR(components::Array{String,1}; userlocations::Array{String,1}=String[], idealmodel = BasicIdeal,verbose=false)
    params = getparams(components, ["properties/critical.csv", "SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"].values
    Tc = params["Tc"]
    Tc_ = Tc.values
    acentricfactor = params["w"]
    idealmodel = idealmodelselector(idealmodel, components,verbose=verbose)
    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 0.457235*R̄^2*Tc_^2/pc/1e6), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. 0.077796*R̄*Tc_/pc/1e6))

    packagedparams = PRParam(a, b, acentricfactor, Tc)
    return PR(packagedparams,idealmodel)
end

function cubic_ab(model::PRModel,T,x)
    a = model.params.a.values
    b = model.params.b.values
    ω = model.params.acentricfactor.values
    Tc = model.params.Tc.values
    α = @. (1+(0.37464+1.54226*ω-0.26992*ω^2)*(1-√(T/Tc)))^2
    āᾱ = sum(a .* .√(α * α') .* (x * x'))
    b̄ = sum(b .* (x * x'))
    return āᾱ ,b̄
end

function cubic_abp(model::PRModel, V, T, x)
    x = z/sum(z)
    n = sum(z)
    v = V/n
    āᾱ ,b̄ = cubic_ab(model,T,x)
    _1 = one(b̄)
    denom = evalpoly(v,(-b̄*b̄,2*b̄,_1))
    p = R̄*T/(v-b̄) - āᾱ /denom
    return a,b,p
end

function cubic_poly(model::PRModel,p,t,x)
    a,b = cubic_ab(model,t,x)
    RT⁻¹ = 1/(R̄*T)
    A = a*p*RT⁻¹*RT⁻¹
    B = b*p*RT⁻¹
    k₀ = B*(B*(B+1.0)-A)
    k₁ = -B*(2*B+1.0) + A
    k₂ = -1.0
    k₃ = 1.0
    return (k₀,k₁,k₂,k₃)
end

function a_resx(model::PRModel, v, T, x)

    a,b = cubic_ab(model,T,x)
    Δ1 = 1+√2
    Δ2 = 1-√2
    ΔPRΔ = 2*√2
    RT⁻¹ = 1/(R̄*T)
    ρ = 1/v
    return -log(1-b*ρ) - a*RT⁻¹*log((Δ1*b*ρ+1)/(Δ2*b*ρ+1))/(ΔPRΔ*b)
    
    
    #return -log(V-n*b̄) + āᾱ/(R̄*T*b̄*2^(3/2)) * log((2*V-2^(3/2)*b̄*n+2*b̄*n)/(2*V+2^(3/2)*b̄*n+2*b̄*n))
end
