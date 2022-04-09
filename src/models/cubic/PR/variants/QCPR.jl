"""
    QCPR(components::Vector{String}; idealmodel=BasicIdeal,
        userlocations=String[], 
        ideal_userlocations=String[],
        alpha_userlocations = String[],
        mixing_userlocations = String[],
        activity_userlocations = String[],
        translation_userlocations = String[],
        verbose=false)

Quantum-corrected Peng Robinson equation of state. it uses the following models:

- Translation Model: [`ConstantTranslation`](@ref)
- Alpha Model: [`TwuAlpha`](@ref)
- Mixing Rule Model: [`QCPRRule`](@ref)

## References

1. Aasen, A., Hammer, M., Lasala, S., Jaubert, J.-N., & Wilhelmsen, Ø. (2020). Accurate quantum-corrected cubic equations of state for helium, neon, hydrogen, deuterium and their mixtures. Fluid Phase Equilibria, 524(112790), 112790. doi:10.1016/j.fluid.2020.112790

"""
function QCPR(components::Vector{String}; idealmodel=BasicIdeal,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

    userlocations = vcat(alpha_userlocations,translation_userlocations,userlocations)

    params = getparams(components, ["cubic/QCPR/QCPR_critical.csv", "cubic/QCPR/QCPR_unlike.csv","cubic/QCPR/Twu_QCPR.csv","cubic/QCPR/QCPR_translation.csv"]; userlocations=userlocations, verbose=verbose)
    k  = params["k"]
    pc = params["pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    init_mixing = init_model(QCPRRule,components,nothing,mixing_userlocations,activity_userlocations,verbose)
    a,b = ab_premixing(PR,init_mixing,Tc,pc,k)
    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)

    M = params["M"]
    N = params["N"]
    L = params["L"]
    packagedparams = TwuAlphaParam(M,N,L)
    init_alpha = TwuAlpha(packagedparams, verbose=verbose)

    c = params["c"]
    packagedparams = ConstantTranslationParam(c)
    init_translation = ConstantTranslation(packagedparams, verbose=verbose)

    icomponents = 1:length(components)
    packagedparams = PRParam(a,b,Tc,pc,Mw)
    references = String["10.1021/I160057A011"]
    model = PR(components,icomponents,init_alpha,init_mixing,init_translation,packagedparams,init_idealmodel,references)
    return model
end

export QCPR


const QCPRModel = PR{T,TwuAlpha,QCPRRule,ConstantTranslation} where T

function cubic_ab(model::QCPRModel,V,T,z=SA[1.0],n=sum(z))
    invn2 = (one(n)/n)^2
    a = model.params.a.values
    b = model.params.b.values
    T = T*float(one(T))
    α = @f(α_function,model.alpha)
    c = @f(translation,model.translation)
    if length(z)>1
    ā,b̄,c̄ = @f(mixing_rule,model.mixing,α,a,b,c)
    else 
        ā = a[1,1]*α[1]
        A = model.mixing.params.A.values[1,1]
        B = model.mixing.params.B.values[1,1]
        Tc = model.params.Tc.values[1]
        β = (1 + A/(T + B))^3 / (1 + A/(Tc + B))^3
        b̄ = b[1,1]*β
        c̄ = c[1]
    end
    return ā ,b̄, c̄
end

function lb_volume(model::QCPRModel,z=SA[1.0])
    A = model.mixing.params.A.values
    B = model.mixing.params.B.values
    l = mixing_model.params.l.values
    Tc = model.params.Tc.values
    n = sum(z)
    invn = (one(n)/n)
    c = model.translation.params.c.values
    c̄ = dot(c,z)
    b̄ = zero(first(z))
    for i in 1:length(z)
        zi = z[i]
        zi2 = zi^2
        Bi = B[i]
        Ai = A[i]
        βi = (1 + Ai/Bi)^3 / (1 + Ai/(Tc[i] + Bi))^3
        bqi = βi*b[i,i]
        b̄ += bqi*zi2
        for j in 1:(i-1)
            zij = zi*z[j]
            Bj = B[j]
            Aj = A[j]
            βj = (1 + Aj/Bj)^3 / (1 + Aj/(Tc[j] + Bj))^3
            bqj = βj*b[j,j]
            b̄ += zij*(bqi+bqj)*(1-l[i,j]) #2 * zij * 0.5(bi + bj)
        end
    end
    return invn*(b̄ - c̄)
end

