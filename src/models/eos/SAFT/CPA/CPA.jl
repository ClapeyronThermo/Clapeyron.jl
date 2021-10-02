struct CPAParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    c1::SingleParam{Float64}
    Tc::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
    Mw::SingleParam{Float64}
end

abstract type CPAModel <: EoSModel end
@newmodel CPA CPAModel CPAParam

export CPA
function CPA(components; idealmodel=BasicIdeal, userlocations=String[], ideal_userlocations=String[], verbose=false)
    params,sites = getparams(components, ["SAFT/CPA", "properties/molarmass.csv","properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    Mw  = params["Mw"]
    k  = params["k"]
    Tc = params["Tc"]
    c1 = params["c1"]
    params["a"].values .*= 1E-1
    params["b"].values .*= 1E-3
    a  = epsilon_LorentzBerthelot(params["a"], k)
    b  = sigma_LorentzBerthelot(params["b"])
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    packagedparams = CPAParam(a, b, c1, Tc, epsilon_assoc, bondvol,Mw)
    references = ["10.1021/ie051305v"]

    model = CPA(packagedparams, sites, idealmodel; ideal_userlocations=ideal_userlocations, references=references, verbose=verbose)
    return model
end

function a_res(model::CPAModel, V, T, z)
    return @f(a_SRK) + @f(a_assoc) + log(V)  # + f(x)
end

function a_SRK(model::CPAModel, V, T, z)
    x = z/∑(z)
    n = ∑(z)
    a = model.params.a.values
    b = model.params.b.values
    Tc = model.params.Tc.values
    c1 = model.params.c1.values

    α = @. (1+c1*(1-√(T/Tc)))^2

    āᾱ = ∑(a .* .√(α * α') .* (x * x'))
    b̄ = ∑(b .* (x * x'))

    return -log(V-n*b̄) - āᾱ/(R̄*T*b̄)*log(1+n*b̄/V)
end
#same as SRK
function ab_consts(model::CPAModel)
    Ωa =  1/(9*(2^(1/3)-1))
    Ωb = (2^(1/3)-1)/3
    return Ωa,Ωb
end

# function cubic_ab(model::CPAModel, T, z)
#     x = z/∑(z)
#     n = ∑(z)
#     a = model.params.a.values
#     b = model.params.b.values
#     Tc = model.params.Tc.values
#     c1 = model.params.c1.values

#     α = @. min((1+c1*(1-√(T/Tc)))^2,one(T))

#     āᾱ = ∑(a .* .√(α * α') .* (x * x'))
#     b̄ = ∑(b .* (x * x'))

#     return āᾱ ,b̄
#end

function Δ(model::CPAModel, V, T, z, i, j, a, b)
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b] * 1e2/R̄
    βijab = model.params.bondvol.values[i,j][a,b] * 1e-3
    Σz = sum(z)
    b = model.params.b.values
    b̄ = dot(z,b,z)
    η = b̄/(4*V*Σz)
    g = (1-η/2)/(1-η)^3
    bij = (b[i,i]+b[j,j])/2
    return g*(exp(ϵ_associjab/T)-1)*βijab*bij/N_A
end