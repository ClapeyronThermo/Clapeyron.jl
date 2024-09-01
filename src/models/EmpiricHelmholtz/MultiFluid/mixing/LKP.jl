
struct LKPMixingParam <: EoSParam
    k::PairParam{Float64}
end

@newmodelsimple LKPMixing LKPMixingParam MixingRule

function vT_scale(model::MultiFluid,mixing::LKPMixing,V,T,z,∑z)
    Vc = model.params.Vc.values
    Tc = model.params.Tc.values
    nc = length(model)
    T̄ = zero(1. + first(z))
    V̄η = zero(T̄)
    V̄ = zero(T̄)
    k = model.params.k.values
    for i in 1:nc
        Vci,Tci,zi = Vc[i],Tc[i],z[i]
        Vciη = sqrt(sqrt(Vci)) #Vci^0.25
        T̄ += zi*zi*Tci*Vciη
        V̄ += zi*zi*Vci
        for j in 1:(i-1)
            Vcj,Tcj,zj = Vc[j],Tc[j],z[j]
            Vcij = 0.125*(cbrt(Vci) + cbrt(Vcj))^3
            Vcijη = sqrt(sqrt(Vcij)) #Vci^0.25
            Tcij = sqrt(Tci*Tcj)*(1 - k[i,j])
            T̄ += 2*zi*zj*Tcij*Vcijη
            V̄ += 2*zi*zj*Vcij
        end
    end
    V̄ = V̄/∑z/∑z
    V̄η = V̄^0.25
    T̄ = T̄/∑z/∑z/V̄η
    return T̄,V̄
end
