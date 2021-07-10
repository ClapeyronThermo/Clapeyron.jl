struct vdWParam <: EoSParam
    Tc::SingleParam{Float64}
    pc::SingleParam{Float64}
    Mw::SingleParam{Float64}
    a::PairParam{Float64}
    b::PairParam{Float64}
end

abstract type vdWModel <: ABCubicModel end
@newmodel vdW vdWModel vdWParam

export vdW
function vdW(components::Array{String,1}; userlocations::Array{String,1}=String[], verbose=false,idealmodel=BasicIdeal)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv","SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    Mw = params["Mw"]
    k = params["k"]
    _pc = params["pc"]
    _Tc = params["Tc"]
    pc = _pc.values
    Tc = _Tc.values

    a = epsilon_LorentzBerthelot(SingleParam(params["pc"], @. 27/64*R̄^2*Tc^2/pc), k)
    b = sigma_LorentzBerthelot(SingleParam(params["pc"], @. 1/8*R̄*Tc/pc))

    packagedparams = vdWParam(_Tc,_pc,Mw,a,b)
    model = vdW(packagedparams,idealmodel)
    return model
end

function ab_consts(::Type{<:vdWModel})
    Ωa =  27/64
    Ωb =  1/8
    return Ωa,Ωb
end

function cubic_ab(model::vdWModel,T,z=SA[1.0],Σn=sum(z))
    inv2Σn = (one(T)/Σn)^2
    a = model.params.a.values
    b = model.params.b.values
    ā = dot(z,Symmetric(a),z)*inv2Σn
    b̄ = dot(z,Symmetric(b),z)*inv2Σn
    return ā,b̄
end

function cubic_abp(model::vdWModel, V, T, z=@SVector [1.0])
    n = ∑(z)
    v = V/n
    a,b = cubic_ab(model,T,z,n)
    p = R̄*T/(v-b) - a/(v*v)
    return a,b,p
end

function cubic_poly(model::vdWModel,p,T,z)
    x = z/sum(z)
    a,b = cubic_ab(model,T,x)
    _1 = one(a+b)
    RT⁻¹ = 1/(R̄*T)
    A = a*p*RT⁻¹*RT⁻¹
    B = b*p*RT⁻¹
    return [-A*B, A, -B-_1, _1]
end



function a_res(model::vdWModel, V, T, z)
    n = sum(z)
    ā,b̄ = cubic_ab(model,T,z,n)
    RT⁻¹ = 1/(R̄*T)
    ρ = n/V
    return -log(1-b̄*ρ) - ā*ρ*RT⁻¹
    # -log(1-n*b/B)-a*n/(R̄*T*B)
    return -log(1-n*b̄/V) - ā*n/(R̄*T*V)
    #return -log(V-n*b̄) - ā*n/(R̄*T*V) + log(V)
end

cubic_zc(::vdWModel) = 3/8
