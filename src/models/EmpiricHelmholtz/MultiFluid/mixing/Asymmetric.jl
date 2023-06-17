struct AsymmetricMixingParam <: EoSParam
    gamma_T::PairParam{Float64}
    gamma_v::PairParam{Float64}
    beta_T::PairParam{Float64}
    beta_v::PairParam{Float64}
end

@newmodelsimple AsymmetricMixing MixingRule AsymmetricMixingParam

function AsymmetricMixing(components;userlocations = String[],verbose = false)
    params = getparams(components,["Empiric/mixing/AsymmetricMixing/AsymmetricMixing_unlike.csv"]; asymmetricparams = ["beta_v","beta_T"],userlocations=userlocations, verbose=verbose)
    beta_v = params["beta_v"]
    gamma_v = params["gamma_v"]
    beta_T = params["beta_T"]
    gamma_T = params["gamma_T"]
    mirror_pair!(beta_T,inv)
    mirror_pair!(beta_v,inv)
    pkgparams = AsymmetricMixingParam(gamma_T,gamma_v,beta_T,beta_v)
    return AsymmetricMixing(pkgparams,verbose = verbose)
end

function recombine_mixing!(model::EmpiricMultiFluid,mixing::AsymmetricMixing)
    Vc = model.params.Vc.values
    Tc = model.params.Tc.values
    nothing
    #TODO
end

function v_scale(model::EmpiricMultiFluid,z,mixing::AsymmetricMixing,∑z)
    vc = model.params.Vc.values
    res = mixing_rule_asymmetric(
        mix_mean3,
        _gerg_asymmetric_mix_rule,
        z,
        vc,
        mixing.params.gamma_v.values,
        mixing.params.beta_v.values,
    )
    return res/(∑z*∑z)
end

function T_scale(model::EmpiricMultiFluid,z,mixing::AsymmetricMixing,∑z)
    Tc = model.params.Tc.values
    #isone(length(z)) && return only(Tc)
    return mixing_rule_asymmetric(
        mix_geomean,
        _gerg_asymmetric_mix_rule,
        z,
        Tc,
        mixing.params.gamma_T.values,
        mixing.params.beta_T.values,
    )/(∑z*∑z)
end

"""
    mixing_rule_asymmetric(op, op_asym, x, p, A, A_asym)

returns an efficient implementation of:
` sum(A[i,j] * x[i] * x[j] * op(p[i],p[j]) * op_asym(x[i],x[j],A_asym[i,j])) for i = 1:n , j = 1:n)`
where `op(p[i],p[j]) == op(p[j],p[i])` , op_asym doesn't follow this symmetry.

""" 
function mixing_rule_asymmetric(op, op_asym, x, p, A, A_asym)
    N = length(x)
    checkbounds(A, N, N)
    checkbounds(A_asym, N, N)
    @boundscheck checkbounds(p, N)
    @inbounds begin
        res1 = zero(eltype(x))
        for i = 1:N
            xi = x[i]
            xi != 0 && begin
                p_i = p[i]
                res1 += p_i * xi^2
                for j = 1:i - 1
                    res1 += 2*xi*x[j]*op(p_i, p[j])*A[i, j]*op_asym(xi, x[j], A_asym[i, j])
                end
            end
        end
    end
    
    return res1
end

_gerg_asymmetric_mix_rule(xi, xj, b) = b * (xi + xj) / (xi * b^2 + xj)

export AsymmetricMixing