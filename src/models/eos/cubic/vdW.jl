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

function cubic_ab(model::vdWModel,T,x)
    a = model.params.a.values
    b = model.params.b.values
    ā = sum(a .* (x * x'))
    b̄ = sum(b .* (x * x'))
    return ā,b̄
end

function cubic_abp(model::vdWModel, V, T, z=@SVector [1.0])
    x = z/sum(z)
    n = sum(z)
    v = V/n
    a,b = cubic_ab(model,T,x)
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



function a_resx(model::vdWModel, v, T, x)
    a,b = cubic_ab(model,T,x)
    RT⁻¹ = 1/(R̄*T)
    ρ = 1/v
    -log(1-b*ρ) - a*ρ*RT⁻¹
    #return -log(V-n*b̄) - ā*n/(R̄*T*V)
end

cubic_zc(::vdWModel) = 3/8
