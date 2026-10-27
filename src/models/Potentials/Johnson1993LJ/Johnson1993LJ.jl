struct Johnson1993LJ <: LJPotentialModel end

function a_potential(model::Johnson1993LJ, V, T̄)
    ρ̄  = 1/V
    γ = 3
    F = exp(-γ*ρ̄*ρ̄)
    x = Johnson1993LJConstants
    T2 = T̄*T̄
    T3 = T2*T̄
    T4 = T2*T2
    T_inv = 1/T̄
    T_inv2 = 1/T2
    T_inv3 = 1/T3
    T_inv4 = 1/T4

    a1 = x[1]*T̄ + x[2]*sqrt(T̄) + x[3] + x[4]*T_inv + x[5]*T_inv2
    a2 = x[6]*T̄ + x[7] + x[8]*T_inv + x[9]*T_inv2
    a3 = x[11] + x[10]*T̄ + x[12]*T_inv
    a4 = x[13]
    a5 = x[14]*T_inv + x[15]*T_inv2
    a6 = x[16]*T_inv
    a7 = x[17]*T_inv + x[18]*T_inv2
    a8 = x[19]*T_inv2
    
    b1 = x[20]*T_inv2 + x[21]*T_inv3
    b2 = x[22]*T_inv2 + x[23]*T_inv4
    b3 = x[24]*T_inv2 + x[25]*T_inv3
    b4 = x[26]*T_inv2 + x[27]*T_inv4
    b5 = x[28]*T_inv2 + x[29]*T_inv3
    b6 = x[30]*T_inv2 + x[31]*T_inv3 + x[32]*T_inv4
    G1 = (1-F)/(2γ)
    ρ̄2 = ρ̄*ρ̄
    ρ̄3 = ρ̄*ρ̄2
    ρ̄4 = ρ̄2*ρ̄2
    ρ̄5 = ρ̄3*ρ̄2
    ρ̄6 = ρ̄3*ρ̄3
    ρ̄8 = ρ̄4*ρ̄4
    ρ̄10 = ρ̄5*ρ̄5
    G2 = -(F*ρ̄2 - 2*G1) / 2γ
    G3 = -(F*ρ̄4 - 4*G2) / 2γ
    G4 = -(F*ρ̄6 - 6*G3) / 2γ
    G5 = -(F*ρ̄8 - 8*G4) / 2γ
    G6 = -(F*ρ̄10 - 10*G5) / 2γ
    bG = b1*G1 + b2*G2 + b3*G3 + b4*G4 + b5*G5 + b6*G6
    ā = (a1,a2/2,a3/3,a4/4,a5/5,a6/6,a7/7,a8/8)
    return (evalpoly(ρ̄,ā)*ρ̄ + bG)*T_inv
end

const Johnson1993LJConstants = (0.8623085097507421, 2.976218765822098,-8.402230115796038,0.1054136629203555,
      -0.8564583828174598, 1.582759470107601,0.7639421948305453, 1.753173414312048,
      2.798291772190376e3,-4.8394220260857657e-2,0.9963265197721935,-3.698000291272493e1,
      2.084012299434647e1, 8.305402124717285e1,-9.574799715203068e2,-1.477746229234994e2,
      6.398607852471505e1, 1.603993673294834e1, 6.805916615864377e1,-2.791293578795945e3,
      -6.245128304568454, -8.116836104958410e3, 1.488735559561229e1,-1.059346754655084e4,
      -1.131607632802822e2,-8.867771540418822e3,-3.986982844450543e1,-4.689270299917261e3,
      2.593535277438717e2,-2.694523589434903e3,-7.218487631550215e2, 1.721802063863269e2)


lb_volume(::Johnson1993LJ,T,z) = sum(z)/0.9 #maximum density was 0.9
crit_pure(::Johnson1993LJ) = (1.3130000571792064, 0.1299353770567644, 3.225806692394338)
has_fast_crit_pure(::Johnson1993LJ) = true

"""
    Johnson1993LJ <: LJPotentialModel
    Johnson1993LJ()

## Description

Johnson equation of state for the Lennard-Jones fluid. Default reference LJ model for [`softSAFT`](@ref)

## References
1. Johnson, J. K., Zollweg, J. A., & Gubbins, K. E. (1993). The Lennard-Jones equation of state revisited. Molecular physics, 78(3), 591–618. [doi:10.1080/00268979300100411](https://doi.org/10.1080/00268979300100411)
"""
Johnson1993LJ

export Johnson1993LJ