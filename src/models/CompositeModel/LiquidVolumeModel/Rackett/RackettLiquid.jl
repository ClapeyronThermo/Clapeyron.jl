abstract type RackettLiquidModel <: LiquidVolumeModel end


struct RackettLiquidParam <: EoSParam
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Zc::SingleParam{Float64}
end

@newmodelsimple RackettLiquid RackettLiquidModel RackettLiquidParam

function RackettLiquid(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = params["acentricfactor"]
    Tc = params["Tc"]
    Pc = params["Pc"]
    vc = params["Vc"]
    _zc = Pc.values .* vc.values ./ (R̄ .* Tc.values)
    Zc = SingleParam("Critical Compressibility factor",components,_zc) 
    packagedparams = RackettLiquidParam(Tc,Pc,Zc)
    model = RackettLiquid(packagedparams)
    return model
end

function volume_impl(model::RackettLiquidModel,p,T,z=SA[1.0],phase=:unknown,threaded=false,vol0 = 0.0)
    tci = model.params.Tc.values
    pci = model.params.Pc.values
    zci = model.params.Zc.values
    Zcm = zero(eltype(z))
    a = zero(eltype(z))
    b = zero(eltype(z)) 
    ∑z = sum(z)
    checkbounds(tci,length(z))
    for i ∈ @comps
        zi = z[i]
        Tcᵢ = tci[i]
        Pcᵢ = pci[i]
        bi = (R̄*Tcᵢ)/Pcᵢ
        ai = (R̄*Tcᵢ)*bi
        zii = zi*zi
        a += zii*ai
        b += zii*bi
        Zcm += z[i]*zci[i]
        for j in 1:(i-1)
            zj = z[j]
            Tcⱼ = tci[i]
            Pcⱼ = pci[i]
            bj = (R̄*Tcⱼ)/Pcⱼ
            aj = (R̄*Tcⱼ)*bi
            aij = sqrt(ai*aj)
            bij = 0.5*(bi+bj)
            zij = zi*zj
            b += 2*zij*bij
            a += 2*zij*aij
        end
    end
    Tcm = a/b/R̄
    Pcm_inv = (b/(∑z*∑z))/(R̄*Tcm)
    Zcm = Zcm/∑z
    Tr = T/Tcm
    return ∑z*R̄*Tcm*Pcm_inv*Zcm^(1+(1-Tr)^(2/7))
end

function volume_impl(model::RackettLiquidModel,p,T,z::SingleComp,phase=:unknown,threaded=false,vol0 = 0.0)
    Tc = only(model.params.Tc.values)
    Pc = only(model.params.Pc.values)
    Pc_inv = 1/Pc
    Zc = only(model.params.Zc.values)
    ∑z = only(z)
    Tr = T/Tc
    return ∑z*R̄*Tc*Pc_inv*Zc^(1+(1-Tr)^(2/7))
end


function YamadaGunnLiquid(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    acentricfactor = params["acentricfactor"]
    Tc = params["Tc"]
    Pc = params["Pc"]
    _zc = 0.29056 .- 0.08775 .* acentricfactor.values
    Zc = SingleParam("Critical Compressibility factor",components,_zc) 
    packagedparams = RackettLiquidParam(Tc,Pc,Zc)
    model = RackettLiquid(packagedparams)
    return model
end
    
export YamadaGunnLiquid,RackettLiquid