struct solidsoftSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
end

abstract type solidsoftSAFTModel <: softSAFTModel end
@newmodel solidsoftSAFT solidsoftSAFTModel solidsoftSAFTParam
default_references(::Type{solidsoftSAFT}) = ["10.1080/002689797170707","10.1080/00268979300100411","10.1080/00268976.2023.2204150"]
default_locations(::Type{solidsoftSAFT}) = ["SAFT/softSAFT/solidsoftSAFT","properties/molarmass.csv"]

function transform_params(::Type{solidsoftSAFT},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

"""
    solidsoftSAFTModel <: SAFTModel

    solidsoftSAFT(components; 
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
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Parameter (no units)
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

Soft SAFT equation of state for the solid phase.

## References
1. Johnson, J. K., Zollweg, J. A., & Gubbins, K. E. (1993). The Lennard-Jones equation of state revisited. Molecular physics, 78(3), 591–618. [doi:10.1080/00268979300100411](https://doi.org/10.1080/00268979300100411)
1. FELIPE J. BLAS and LOURDES F. VEGA. (1997). Thermodynamic behaviour of homonuclear and heteronuclear Lennard-Jones chains with association sites from simulation and theory. Molecular physics, 92(1), 135–150. [doi:10.1080/002689797170707](https://doi.org/10.1080/002689797170707)
3. Ramírez-Carpio, V., Galindo, A., & Gil-Villegas, A. (2023). Modelling the solid–liquid–vapour phase behaviour of n -alkanes in a TPT-1 framework. Molecular Physics, 121(19–20). [doi:10.1080/00268976.2023.2204150](https://doi.org/10.1080/00268976.2023.2204150)
"""
solidsoftSAFT

export solidsoftSAFT

recombine_impl!(model::solidsoftSAFTModel) = recombine_saft!(model)
is_solid(model::solidsoftSAFTModel) = true

function x0_volume_solid(model::solidsoftSAFTModel,T,z)
    v_lb = lb_volume(model,z)
    return v_lb*1.2
end

function data(model::solidsoftSAFTModel,V,T,z)
    σ3,ϵ̄,m̄ = @f(σϵ_m_vdw1f)
    ∑z = sum(z)
    N = N_A*∑z
    ρS = N/V*m̄
    ρ̄  = ρS*σ3
    return σ3,ϵ̄,m̄,ρ̄ 
end

function a_res(model::solidsoftSAFTModel, V, T, z)    
    _data = @f(data)
    return @f(a_LJ,_data) + @f(a_chain,_data)
end

function a_LJ(model::solidsoftSAFTModel, V, T, z,_data = @f(data))
    σ3,ϵ̄,m̄,ρ̄  = _data
    T̄ = T/ϵ̄
    a = solidsoftSAFTconsts.a_ref
    b = solidsoftSAFTconsts.b
    ustat = -14.45392093*ρ̄^2+6.065940096*ρ̄^4
    Uah = -sum(a[n+1,m-1]*ρ̄^n*T̄^(m-1)/(m-1) for n ∈ 0:2 for m ∈ 2:5)
    B = sum(b[n]/n*ρ̄^n for n ∈ 1:4)
    return m̄*(-23.3450759+ustat/T̄+Uah-3/2*log(T̄)+B)
end

function a_chain(model::solidsoftSAFTModel, V, T, z,_data = @f(data))
    σ3,ϵ̄,m̄,ρ̄  = _data
    return -log(@f(g_LJ,_data)*exp(-ϵ̄/T))*(m̄-1)
end

function g_LJ(model::solidsoftSAFTModel, V, T, z ,_data = @f(data))
    σ3,ϵ̄,m̄,ρ̄  = _data
    T̄ = T/ϵ̄
    a = solidsoftSAFTconsts.a
    gLJ = 1+sum(a[i,j]*ρ̄^i*T̄^(1-j) for i ∈ 1:5 for j ∈ 1:5)
    return gLJ
end


const solidsoftSAFTconsts =
(
 a_ref = [ -8.2151768  12.070686 -6.6594615  1.3211582;
            13.404069 -20.632066  11.564825 -2.3064801;
           -5.5481261  8.8465978 -5.0258631  1.0070066 ],

 a = [-11.632  37.706 -140.655   52.675   1.019;
       86.742 -40.865  335.679 -108.881  -17.970;
     -131.434 -190.01 -110.953   -2.908   48.886;
       69.219 311.947 -197.314  114.210  -47.051;
      -10.560 -120.436  112.935  -54.753   15.058],

 b = [69.833875,-132.86963, 97.438593,-25.848057],
)
