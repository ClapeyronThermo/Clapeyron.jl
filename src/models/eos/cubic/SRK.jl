struct SRKParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    acentricfactor::SingleParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type SRKModel <: ABCubicModel end
@newmodel SRK SRKModel SRKParam

export SRK
function SRK(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false,idealmodel=BasicIdeal)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    k  = params["k"]
    _pc = params["pc"]
    pc = _pc.values


    Tc = params["Tc"]
    Tc_ = Tc.values
    acentricfactor = params["w"]
    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 1/(9*(2^(1/3)-1))*R̄^2*Tc_^2/pc), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. (2^(1/3)-1)/3*R̄*Tc_/pc))

    packagedparams = SRKParam(a, b, acentricfactor, Tc,_pc,Mw)
    model = SRK(packagedparams,idealmodel)
    return model
end

function ab_consts(::Type{<:SRKModel})
    Ωa =  1/(9*(2^(1/3)-1))
    Ωb = (2^(1/3)-1)/3
    return Ωa,Ωb
end

function cubic_ab(model::SRKModel,T,x)
    a = model.params.a.values
    b = model.params.b.values
    ω = model.params.acentricfactor.values
    Tc = model.params.Tc.values
    _1 = one(T+first(x))
    #root of α
    αx = @. sqrt(min((_1+(0.480+1.547*ω-0.176*ω^2)*(1-√(T/Tc)))^2,_1)) * x
    #αx = (αi*xi for (αi,xi) in zip(α,x)) #lazy iterator
    #αx .*= x
    #āᾱ = sum(a .* .√(α * α') .* (x * x'))
    āᾱ =dot(αx, Symmetric(a), αx)
    b̄ = dot(x, Symmetric(b), x)
    return āᾱ ,b̄
end

function cubic_abp(model::SRKModel, V, T, x)
    a,b = cubic_ab(model,T,x)
    p =  R̄*T/(V-b) - a/((V+b)*V)
    return a,b,p
end

function cubic_poly(model::SRKModel,p,T,z)
    x = z/sum(z)
    a,b = cubic_ab(model,T,x)
    RT⁻¹ = 1/(R̄*T)
    A = a*p* RT⁻¹* RT⁻¹
    B = b*p* RT⁻¹
    _1 = one(a)
    return [-A*B, -B*(B+_1) + A, -_1, _1]
end

function a_resx(model::SRKModel, v, T, x)

    ā,b̄ = cubic_ab(model,T,x)
    ρ = 1/v
    RT⁻¹ = 1/(R̄*T)
    return -log(1-b̄*ρ) - ā*RT⁻¹*log(b̄*ρ+1)/b̄
    #=
    x = z/sum(z)
    n = sum(z)
    a = model.params.a.values
    b = model.params.b.values
    ω = model.params.acentricfactor.values
    Tc = model.params.Tc.values

    α = @. (1+(0.480+1.547*ω-0.176*ω^2)*(1-√(T/Tc)))^2

    āᾱ = sum(a .* .√(α * α') .* (x * x'))
    b̄ = sum(b .* (x * x'))
    return -log(V-n*b̄) - āᾱ/(R̄*T*b̄)*log(1+n*b̄/V)
    =#
end

cubic_zc(::SRKModel) = 1/3
