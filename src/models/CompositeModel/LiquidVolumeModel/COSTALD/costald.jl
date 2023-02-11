abstract type COSTALDModel <: LiquidVolumeModel end


struct COSTALDParam <: EoSParam
    Tc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    acentricfactor::SingleParam{Float64}
end

@newmodelsimple COSTALD COSTALDModel COSTALDParam

function COSTALD(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components, ["properties/critical.csv"]; userlocations=userlocations, verbose=verbose)
    Tc = params["Tc"]
    Vc = params["Vc"]
    acentricfactor = params["acentricfactor"]
    packagedparams = COSTALDParam(Tc,Vc,acentricfactor)
    model = COSTALD(packagedparams;verbose)
    return model
end

function volume_impl(model::COSTALDModel,p,T,z=SA[1.0],phase=:unknown,threaded=false,vol0 = nothing)
    Tci = model.params.Tc.values
    Vci = model.params.Vc.values
    ωi  = model.params.acentricfactor.values
        
    Vc1 = zero(eltype(z))
    Vc23 = zero(eltype(z))
    Vc13 = zero(eltype(z))
    ω = zero(eltype(z))  
    Tc = zero(eltype(z))
    checkbounds(Tci,length(z))
    for i ∈ @comps
        zi = z[i]
        ω += zi*ωi[i]
        Vcii = Vci[i]
        Vcii13 = cbrt(Vcii)
        Vcii23 = Vcii13*Vcii13
        Vc1 += zi*Vcii
        Vc13 += zi*Vcii13
        Vc23 += zi*Vcii23
        VTi = Vcii*Tci[i]
        Tc += zi*zi*VTi
        for j ∈ 1:i-1
            VTj = Vci[j]*Tci[j]
            Tc += 2*zi*z[j]*sqrt(VTi*VTj)
        end
    end

    Vc = 0.25*(Vc1 + 3*Vc13*Vc23)
    Tc = Tc/Vc

    Tr = T/Tc
    τ = 1.0 - Tr
    τcbrt = cbrt(τ)
    Vδ = evalpoly(Tr,(-0.296123,0.386914,-0.0427258,-0.0480645))/(Tr - 1.00001)
    V0 = evalpoly(τcbrt,(1.0, -1.52816,1.43907,-0.81446,0.190454))
    return Vc*V0*(1.0 - ω*Vδ)
end

function volume_impl(model::COSTALDModel,p,T,z::SingleComp,phase=:unknown,threaded=false,vol0 = nothing)
    Tc = model.params.Tc.values |> only
    Vc = model.params.Vc.values |> only
    ω  = model.params.acentricfactor.values |> only

    Tr = T/Tc
    τ = 1.0 - Tr
    τcbrt = cbrt(τ)
    Vδ = evalpoly(Tr,(-0.296123,0.386914,-0.0427258,-0.0480645))/(Tr - 1.00001)
    V0 = evalpoly(τcbrt,(1.0, -1.52816,1.43907,-0.81446,0.190454))
    return Vc*V0*(1.0 - ω*Vδ)
end

export COSTALD