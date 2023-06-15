struct LorentzBerthelotMixingParam <: MixingRule
    k::PairParam{Float64}
    l::PairParam{Float64}
end

@newmodelsimple LorentzBerthelotMixing MixingRule LorentzBerthelotMixingParam

function LorentzBerthelotMixing(components;userlocations = String[],verbose = false)
    params = getparams(components,["Empiric/mixing/LorentzBerthelotMixing/LorentzBerthelotMixing_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    l = params["l"]
    k = params["k"]
    pkgparams = LorentzBerthelotMixingParam(k,l)
    return LorentzBerthelotMixing(pkgparams, verbose = verbose)
end

function v_scale(model::EmpiricMultiFluid,z,mixing::LorentzBerthelotMixing,∑z)
    Vc = model.params.Tc.values
    l = mixing.params.l
    res = zero(∑z)*1.0
    for i in @comps
        Vci = Vc[i]
        zi = z[i]
        res += zi*zi*Vci
        for j in 1:(i-1)
            Vcj = Vc[j]
            zi = z[j]
            res += zi*zj*mix_mean3(Vci,Vcj,l[i,j])
        end
    end
    return res/(∑z*∑z)
end

function T_scale(model::EmpiricMultiFluid,z,mixing::LorentzBerthelotMixing,∑z)
    Tc = model.params.Tc.values
    k = mixing.params.k
    res = zero(∑z)*1.0
    for i in @comps
        Tci = Vc[i]
        zi = z[i]
        res += zi*zi*Tci
        for j in 1:(i-1)
            Tcj = Vc[j]
            zi = z[j]
            res += zi*zj*mix_geomean(Tci,Tcj,k[i,j])
        end
    end
    return res/(∑z*∑z)
end

export LorentzBerthelotMixing