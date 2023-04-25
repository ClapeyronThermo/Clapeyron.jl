struct PPCSAFTParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    dipole::SingleParam{Float64}
    epsilon_assoc::AssocParam{Float64}
    bondvol::AssocParam{Float64}
end

abstract type PPCSAFTModel <: PCSAFTModel end
@newmodel PPCSAFT PPCSAFTModel PPCSAFTParam

"""
    PPCSAFTModel <: SAFTModel
    PPCSAFT(components; 
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`
## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `dipole`: Single Parameter (`Float64`) - Dipole moment `[D]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume
## Input models
- `idealmodel`: Ideal Model
## Description
Polar Perturbed-Chain SAFT (PPC-SAFT)
## References
1. Gross, J., & Vrabec, J. (2005). An equation-of-state contribution for polar components: Dipolar molecules. AIChE Journal, 52(3), 856-1282. [doi:10.1002/aic.10683](https://doi.org/10.1002/aic.10683)
"""
PPCSAFT

export PPCSAFT
function PPCSAFT(components;
    idealmodel=BasicIdeal,
    userlocations=String[],
    ideal_userlocations=String[],
    verbose=false,
    assoc_options = AssocOptions())
    params,sites = getparams(components, ["SAFT/PCSAFT/PPCSAFT/","properties/molarmass.csv"]; userlocations=userlocations, verbose=verbose)
    segment = params["segment"]
    k = get(params,"k",nothing)
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    bondvol,epsilon_assoc = assoc_mix(bondvol,epsilon_assoc,sigma,assoc_options) #combining rules for association

    dipole = params["dipole"]

    packagedparams = PPCSAFTParam(Mw,segment,sigma,epsilon,dipole,epsilon_assoc,bondvol)
    references = ["10.1021/ie0003887", "10.1021/ie010954d"]

    model = PPCSAFT(packagedparams, sites, idealmodel; ideal_userlocations, references, verbose, assoc_options)
    return model
end

recombine_impl!(model ::PPCSAFTModel) = recombine_saft!(model)

function a_res(model ::PPCSAFTModel, V, T, z)
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data) + @f(a_polar,_data)
end

function a_polar(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    a₂ = @f(a_2,_data)
    a₃ = @f(a_3,_data)
    if a₂ == 0.0
        return 0.0
    else
        return a₂^2/(a₂-a₃)
    end
end

function a_2(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    ∑z = sum(z)
    ρ = N_A*∑z/V
    x = z./∑z
    _a_2 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    μ = model.params.dipole.values
    μ̄² = μ.^2 ./m./k_B*1e-36*(1e-10*1e-3)
    @inbounds for i ∈ @comps
        _J2_ii = @f(J2,i,i,_data)
        _a_2 += x[i]^2*μ̄²[i]^2/σ[i,i]^3/T^2*_J2_ii
        for j ∈ i+1:length(model)
            _J2_ij = @f(J2,i,j,_data)
            _a_2 += 2*x[i]*x[j]*μ̄²[i]*μ̄²[j]/σ[i,j]^3/T^2*_J2_ij
        end
    end
    _a_2 *= -π*ρ
    return _a_2
end

function a_3(model ::PPCSAFTModel, V, T, z, _data=@f(data))
    ∑z = sum(z)
    ρ = N_A*∑z/V
    x = z./∑z
    _a_3 = zero(T+V+first(z))
    m = model.params.segment.values
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    μ = model.params.dipole.values
    μ̄² = μ.^2 ./m./k_B*1e-36*(1e-10*1e-3)
    @inbounds for i ∈ @comps
        _J3_iii = @f(J3,i,i,i,_data)
        _a_3 += x[i]^3*μ̄²[i]^3/σ[i,i]^3/T^3*_J3_iii
        for j ∈ i+1:length(model)
            _J3_iij = @f(J3,i,i,j,_data)
            _a_3 += 3*x[i]^2*x[j]*μ̄²[i]^2*μ̄²[j]/σ[i,j]^2/σ[i,i]/T^3*_J3_iij
            _J3_ijj = @f(J3,i,j,j,_data)
            _a_3 += 3*x[i]*x[j]^2*μ̄²[i]*μ̄²[j]^2/σ[i,j]^2/σ[j,j]/T^3*_J3_ijj
            for k in j+1:length(model)
                _J3_ijk = @f(J3,i,j,k,_data)
                _a_3 += 6*x[i]*x[j]*x[k]*μ̄²[i]*μ̄²[j]*μ̄²[k]/σ[i,j]/σ[i,k]/σ[j,k]/T^3*_J3_ijk
            end
        end
    end
    _a_3 *= -4π^2/3*ρ^2
    return _a_3
end

function J2(model, V, T, z, i, j, _data=@f(data))
    _,_,_,_,η,_ = _data

    corr_a = PPCSAFTconsts.corr_a
    corr_b = PPCSAFTconsts.corr_b

    m = model.params.segment.values
    m̄ = (m[i]*m[j])^(1/2)
    if m̄>2
        m̄=2
    end
    ϵ = model.params.epsilon.values
    res = zero(T+V+first(z))
    @inbounds for k ∈ 1:5
        ii = k-1 
        a1,a2,a3 = corr_a[k]
        aij = a1 + (m̄-1)/m̄*a2 + (m̄-1)/m̄*(m̄-2)/m̄*a3

        b1,b2,b3 = corr_b[k]
        bij = b1 + (m̄-1)/m̄*b2 + (m̄-1)/m̄*(m̄-2)/m̄*b3

        res +=(aij+bij*ϵ[i,j]/T)*η^ii
    end
    return res
end

function J3(model, V, T, z, i, j, k, _data=@f(data))
    _,_,_,_,η,_ = _data

    corr_c = PPCSAFTconsts.corr_c

    m = model.params.segment.values
    m̄ = (m[i]*m[j]*m[k])^(1/3)
    if m̄>2
        m̄=2
    end
    res = zero(η)
    @inbounds for k ∈ 1:5
        ii = k-1 
        c1,c2,c3 = corr_c[k]
        cijk = c1 + (m̄-1)/m̄*c2 + (m̄-1)/m̄*(m̄-2)/m̄*c3

        res +=cijk*η^ii
    end
    return res
end

const PPCSAFTconsts = (
    corr_a =
    [(0.3043504,0.9534641,-1.161008),
    (-0.1358588,-1.8396383,4.5258607),
    (1.4493329,2.013118,0.9751222),
    (0.3556977,-7.3724958,-12.281038),
    (-2.0653308,8.2374135,5.9397575)],

    corr_b =
    [(0.2187939,-0.5873164,3.4869576),
    (-1.1896431,1.2489132,-14.915974),
    (1.1626889,-0.508528,15.372022),
    (0,0,0),
    (0,0,0)],

    corr_c =
    [(-0.0646774,-0.9520876,-0.6260979),
    (0.1975882,2.9924258,1.2924686),
    (-0.8087562,-2.3802636,1.6542783),
    (0.6902849,-0.2701261,-3.4396744),
    (0,0,0)]
)