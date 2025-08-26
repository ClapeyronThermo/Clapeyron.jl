@newmodelsingleton JagerSpanSolidCO2 GibbsBasedModel

"""
    JagerSpanSolidCO2 <: GibbsBasedModel
    JagerSpanSolidCO2()

## Input parameters

None

## Description

Jäger and Span gibbs model for solid carbon dioxide.

## References

1. Jäger, A., & Span, R. (2012). Equation of state for solid carbon dioxide based on the Gibbs free energy. Journal of Chemical and Engineering Data, 57(2), 590–597. [doi:10.1021/je2011677]((https://doi.org/10.1021/je2011677)
"""
JagerSpanSolidCO2

default_references(::Type{JagerSpanSolidCO2}) = ["10.1021/je2011677"]

function eos_g(model::JagerSpanSolidCO2,p,T,z)
    g =
    ϑ,_π = T/150.0,p/101325.0
    Δϑ,Δπ = ϑ - 1,_π - 1
    ϑ2 = ϑ*ϑ
    g₀,g₁,g₂,g₃,g₄,g₅,g₆,g₇,g₈,g₉,g₁₀ = JagerSpanSolidCO2Consts.g
    n = 7
    n̄ = 1 - 1/n
    ĝ₁ = evalpoly(Δϑ,(g₀,g₁,g₂))
    ĝ₂ = g₃*(log((ϑ2 + g₄*g₄)/(1 + g₄*g₄)) - (2*ϑ/g₄)*(atan(ϑ/g₄) - atan(1/g₄)))
    ĝ₃ = g₅*(log((ϑ2 + g₆*g₆)/(1 + g₆*g₆)) - (2*ϑ/g₆)*(atan(ϑ/g₆) - atan(1/g₆)))

    fα = CO2_fα(model,ϑ)
    K = CO2_K(model,ϑ)
    ĝ₄ = g₇*Δπ*(exp(fα) + K*g₈) + g₉*K*((_π + g₁₀)^n̄ - (1 + g₁₀)^n̄)

    ĝ = ĝ₁ + ĝ₂ + ĝ₃ + ĝ₄
    g = ĝ*8.314472*150.0
    return g*sum(z)
end

function CO2_fα(model::JagerSpanSolidCO2,ϑ)
    ϑ2 = ϑ*ϑ
    gα₀,gα₁,gα₂,gα₃,gα₄,gα₅,gα₆,gα₇,gα₈ = JagerSpanSolidCO2Consts.gα
    f1 = gα₀*(ϑ2 - 1)
    f2 = gα₁*log((ϑ2 - gα₂*ϑ + gα₃)/(1 - gα₂ + gα₃))
    f3 = gα₄*log((ϑ2 + gα₂*ϑ + gα₃)/(1 + gα₂ + gα₃))
    f4 = gα₅*(atan((ϑ - gα₆)/gα₇) - atan((1 - gα₆)/gα₇))
    f5 = gα₈*(atan((ϑ + gα₆)/gα₇) - atan((1 + gα₆)/gα₇))
    return f1+f2+f3+f4+f5
end

function CO2_K(model::JagerSpanSolidCO2,ϑ)
    gk₀,gk₁,gk₂ = JagerSpanSolidCO2Consts.gk
    return evalpoly(ϑ,(gk₂,gk₁,gk₀))
end

function gibbsmodel_reference_state_consts(model::JagerSpanSolidCO2)
    return :dH,517950.0,216.592,8875.0
end

function gibbsmodel_reference_state_consts(ice::JagerSpanSolidCO2,water::EmpiricHelmholtzModel)
    return :zero,0.0,0.0,0.0
end

p_scale(model::JagerSpanSolidCO2,z) = 101325.0
T_scale(model::JagerSpanSolidCO2,z) = 150.0

const JagerSpanSolidCO2Consts = (
g = (-2.6385478, 4.5088732, -2.0109135, -2.7976237, 0.26427834, 3.8259935, 0.31711996, 0.0022087195, -1.1289668, 0.0092923982, 3391.4617),
gα = (0.039993365, 0.0023945101, 0.32839467, 0.057918471, 0.0023945101, -0.0026531689, 0.16419734, 0.17594802, 0.0026531689),
gk = (2.2690751e-1,-7.5019750e-2,2.6442913e-1),
)

export JagerSpanSolidCO2
