struct softSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type softSAFTModel <: SAFTModel end
@newmodel softSAFT softSAFTModel softSAFTParam
default_references(::Type{softSAFT}) = ["10.1080/002689797170707","10.1080/00268979300100411"]
default_locations(::Type{softSAFT}) = ["SAFT/softSAFT","properties/molarmass.csv"]
function transform_params(::Type{softSAFT},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

function get_k(model::softSAFTModel)   
    return get_k_geomean(model.params.epsilon)
end

function get_l(model::softSAFTModel)   
    return get_k_mean(model.params.sigma)
end

"""
    softSAFTModel <: SAFTModel

    softSAFT(components; 
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m³]`

## Input models
- `idealmodel`: Ideal Model

## Description

Soft SAFT, with Lennard-Jones function from Johnson et al. (1993)

## References
1. Johnson, J. K., Zollweg, J. A., & Gubbins, K. E. (1993). The Lennard-Jones equation of state revisited. Molecular physics, 78(3), 591–618. [doi:10.1080/00268979300100411](https://doi.org/10.1080/00268979300100411)
2. FELIPE J. BLAS and LOURDES F. VEGA. (1997). Thermodynamic behaviour of homonuclear and heteronuclear Lennard-Jones chains with association sites from simulation and theory. Molecular physics, 92(1), 135–150. [doi:10.1080/002689797170707](https://doi.org/10.1080/002689797170707)
"""
softSAFT

export softSAFT

recombine_impl!(model::softSAFTModel) = recombine_saft!(model)

function lb_volume(model::softSAFTModel,z)
    σ3,ϵ̄,m̄ = σϵ_m_vdw1f(model,1.0,1.0,z)
    return m̄*N_A*σ3*π/6
end


function x0_volume_liquid(model::softSAFTModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.8
end


function data(model::softSAFTModel,V,T,z)
    σ3,ϵ̄,m̄ = σϵ_m_vdw1f(model,V,T,z)
    ∑z = sum(z)
    N = N_A*∑z
    ρS = N/V*m̄
    ρ̄  = ρS*σ3
    return σ3,ϵ̄,m̄,ρ̄ 
end

function a_res(model::softSAFTModel, V, T, z)    
    _data = @f(data)
    return @f(a_LJ,_data) + @f(a_chain) + @f(a_assoc)
end

function a_LJ(model::softSAFTModel, V, T, z,_data = @f(data))
    σ3,ϵ̄,m̄,ρ̄  = _data
    T̄ = T/ϵ̄
    γ = 3
    F = exp(-γ*ρ̄*ρ̄)
    x = softSAFTconsts.x
    T2 = T̄*T̄
    T3 = T2*T̄
    T4 = T2*T2
    T_inv = ϵ̄/T
    T_inv2 = 1/T2
    T_inv3 = 1/T3
    T_inv4 = 1/T4

    a1 = x[1]*T̄ + x[2]*√(T̄) + x[3] + x[4]*T_inv + x[5]*T_inv2
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
    return m̄*(evalpoly(ρ̄,ā)*ρ̄ + bG)*T_inv
end

function a_chain(model::softSAFTModel, V, T, z,_data = @f(data))
    σ3,ϵ̄,m̄,ρ̄  = _data
    return -log(@f(y_LJ,_data))*(m̄-1)
end

function y_LJ(model::softSAFTModel, V, T, z, _data = @f(data))
    σ3,ϵ̄,m̄,ρ̄  = _data
    gLJ = @f(g_LJ,_data)
    return gLJ*exp(-ϵ̄/T)
end

function g_LJ(model::softSAFTModel, V, T, z ,_data = @f(data))
    σ3,ϵ̄,m̄,ρ̄  = _data
    T̄ = T/ϵ̄
    a = softSAFTconsts.a
    gLJ = 1+sum(a[i,j]*ρ̄^i*T̄^(1-j) for i ∈ 1:5 for j ∈ 1:5)
end

function Δ(model::softSAFTModel, V, T, z, i, j, a, b,_data = @f(data))
    σ3,ϵ̄,m̄,ρ̄   = _data
    T̄ = T/ϵ̄
    ϵ_assoc = model.params.epsilon_assoc.values
    κ = model.params.bondvol.values
    b_ = softSAFTconsts.b

    I = sum(b_[i+1,j+1]*ρ̄^i*T̄^j for i ∈ 0:4 for j ∈ 0:4)/3.84/1e4
    return 4π*(exp(ϵ_assoc[i,j][a,b]/T)-1)*κ[i,j][a,b]*I
end


const softSAFTconsts =
(
 x = [ 0.8623085097507421, 2.976218765822098,-8.402230115796038,0.1054136629203555,
      -0.8564583828174598, 1.582759470107601,0.7639421948305453, 1.753173414312048,
      2.798291772190376e3,-4.8394220260857657e-2,0.9963265197721935,-3.698000291272493e1,
      2.084012299434647e1, 8.305402124717285e1,-9.574799715203068e2,-1.477746229234994e2,
      6.398607852471505e1, 1.603993673294834e1, 6.805916615864377e1,-2.791293578795945e3,
      -6.245128304568454, -8.116836104958410e3, 1.488735559561229e1,-1.059346754655084e4,
      -1.131607632802822e2,-8.867771540418822e3,-3.986982844450543e1,-4.689270299917261e3,
      2.593535277438717e2,-2.694523589434903e3,-7.218487631550215e2, 1.721802063863269e2],

 a = [0.49304346593882 2.1528349894745 -15.955682329017 24.035999666294 -8.6437958513990;
      -0.47031983115362 1.1471647487376  37.889828024211 -84.667121491179 39.643914108411;
      5.0325486243620 -25.915399226419 -18.862251310090 107.63707381726 -66.602649735720;
      -7.3633150434385  51.553565337453 -40.519369256098 -38.796692647218 44.605139198378;
      2.9043607296043 -24.478812869291  31.500186765040 -5.3368920371407 -9.5183440180133],

 b = [-0.03915181 0.08450471 0.06889053 -0.01034279  0.5728662e-3;
      -0.5915018  0.9838141 -0.4862279   0.1029708 -0.6919154e-2;
      1.908368  -3.415721   2.124052  -0.4298159  0.02798384;
      -0.7957312  0.7187330 -0.9678804   0.2431675 -0.01644710;
      -0.9399577   2.314054 -0.4877045  0.03932058 -0.1600850e-2],
)

#=
function ϵ_m(model::softSAFTModel, V, T, z)
    comps = @comps
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    m = model.params.segment.values
    return sum(m[i]*m[j]*z[i]*z[j]*σ[i,j]^3*ϵ[i,j] for i ∈ comps for j ∈ comps)/sum(m[i]*m[j]*z[i]*z[j]*σ[i,j]^3 for i ∈ comps for j ∈ comps)
end

function σ_m(model::softSAFTModel, V, T, z)
    comps = @comps
    σ = model.params.sigma.values
    m = model.params.segment.values
    return (sum(m[i]*m[j]*z[i]*z[j]*σ[i,j]^3 for i ∈ comps for j ∈ comps)/sum(m[i]*m[j]*z[i]*z[j] for i ∈ comps for j ∈ comps))^(1/3)
end

function ρ_S(model::softSAFTModel, V, T, z)
    ∑z = ∑(z)
    N = N_A*∑z
    m = model.params.segment.values
    m̄ = dot(z,m)/∑z
    return N/V*m̄
end =#
