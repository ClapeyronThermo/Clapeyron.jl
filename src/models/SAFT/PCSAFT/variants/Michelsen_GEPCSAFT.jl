
abstract type Michelsen_GEPCSAFTModel <: SAFTModel end

const Michelsen_GEPCSAFTParam = PCSAFTParam

struct Michelsen_GEPCSAFTi{T <: IdealModel,γ} <: Michelsen_GEPCSAFTModel
    components::Array{String,1}
    sites::SiteParam
    activity::γ
    params::Michelsen_GEPCSAFTParam
    idealmodel::T
    assoc_options::AssocOptions
    references::Array{String,1}
end

"""
    Michelsen_GEPCSAFT <: SAFTModel

    Michelsen_GEPCSAFT(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    verbose = false,
    assoc_options = AssocOptions())

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter [`A°`]
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy  `[K]`
- `k`: Pair Parameter (`Float64`) (optional) - Binary Interaction Paramater (no units)
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume `[m^3]`

## Model Parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`
- `epsilon_assoc`: Association Parameter (`Float64`) - Reduced association energy `[K]`
- `bondvol`: Association Parameter (`Float64`) - Association Volume

## Input models
- `idealmodel`: Ideal Model
- `activity`: Activity model

## Description

Perturbed-Chain SAFT (PC-SAFT), with Gᴱ mixing rule (using the Michelsen limit).

## References
"""
Michelsen_GEPCSAFT

export Michelsen_GEPCSAFT
function Michelsen_GEPCSAFT(components;
    idealmodel = BasicIdeal,
    activity = UNIFAC,
    userlocations = String[],
    ideal_userlocations = String[],
    activity_userlocations = String[],
    assoc_options = AssocOptions(),
    reference_state = nothing,
    verbose = false)

    
    params = getparams(components, ["SAFT/PCSAFT/PCSAFT_like.csv","SAFT/PCSAFT/PCSAFT_unlike.csv","SAFT/PCSAFT/PCSAFT_assoc.csv"]; userlocations = userlocations, verbose = verbose)

    sites = params["sites"]
    segment = params["segment"]
    k = get(params,"k",nothing)
    Mw = params["Mw"]
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"])
    epsilon = epsilon_LorentzBerthelot(params["epsilon"], k)
    epsilon_assoc = params["epsilon_assoc"]
    bondvol = params["bondvol"]
    packagedparams = Michelsen_GEPCSAFTParam(Mw, segment, sigma, epsilon, epsilon_assoc, bondvol)

    init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
    init_activity = init_model(activity,components,activity_userlocations,verbose)
    references = String["10.1021/acs.iecr.2c03464"]
    
    model = Michelsen_GEPCSAFTi(format_components(components),sites,init_activity,packagedparams,init_idealmodel,assoc_options,references)
    set_reference_state!(model,reference_state;verbose)
    return model 
end

function a_res(model::Michelsen_GEPCSAFTModel, V, T, z)
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function data(model::Michelsen_GEPCSAFTModel,V,T,z)
    _d = @f(d)
    ζ0,ζ1,ζ2,ζ3 = @f(ζ0123,_d)
    m = model.params.segment.values
    m̄ = dot(z, m)/sum(z)
    return (_d,ζ0,ζ1,ζ2,ζ3,m̄)
end

function a_hc(model::Michelsen_GEPCSAFTModel, V, T, z , _data = @f(data))
    _,_,_,_,η,m̄ = _data
    g_hs = (1-η/2)/(1-η)^3
    a_hs = (4η-3η^2)/(1-η)^2
    return m̄*a_hs - (m̄-1)*log(g_hs)
end

function g_hs(model::Michelsen_GEPCSAFTModel, V, T, z,_data = @f(data))
    _,_,_,_,η,_ = _data
    return (1-η/2)/(1-η)^3
end

function a_hs(model::Michelsen_GEPCSAFTModel, V, T, z)
    _,_,_,_,η,_ = _data
    return (4η-3η^2)/(1-η)^2
end

function Δ(model::Michelsen_GEPCSAFTModel, V, T, z, i, j, a, b,_data = @f(data))
    ϵ_associjab = model.params.epsilon_assoc.values[i,j][a,b]
    κijab = model.params.bondvol.values[i,j][a,b]
    σij = model.params.sigma.values[i,j]
    g_hs_ = @f(g_hs,_data)
    return g_hs_*σij^3*(exp(ϵ_associjab/T)-1)*κijab
end

function a_disp(model::Michelsen_GEPCSAFTModel, V, T, z,_data=@f(data))
    di,ζ0,ζ1,ζ2,ζ3,m̄ = _data
    Σz = sum(z)
    m2ϵσ3₁,m2ϵσ3₂ = @f(m2ϵσ3,_data)
    return -2*π*N_A*Σz/V*@f(I,1,_data)*m2ϵσ3₁ - π*m̄*N_A*Σz/V*@f(C1,_data)*@f(I,2,_data)*m2ϵσ3₂
end

function d(model::Michelsen_GEPCSAFTModel, V, T, z)
    ϵᵢᵢ = diagvalues(model.params.epsilon)
    σᵢᵢ = diagvalues(model.params.sigma)
    return σᵢᵢ .* (1 .- 0.12 .* exp.(-3ϵᵢᵢ ./ T))
end


function C1(model::Michelsen_GEPCSAFTModel, V, T, z,_data=@f(data))
    _,_,_,_,η,m̄ = _data
    return (1 + m̄*(8η-2η^2)/(1-η)^4 + (1-m̄)*(20η-27η^2+12η^3-2η^4)/((1-η)*(2-η))^2)^-1
end

function m2ϵσ3(model::Michelsen_GEPCSAFTModel, V, T, z, _data=@f(data))
    Tnum = promote_type(eltype(z), typeof(V), typeof(T))

    N = length(z)
    Nbin = div(N*(N-1),2)

    m_ij = zeros(Tnum, N, N)
    σ_ij = zeros(Tnum, N, N)
    ϵ_ii = zeros(Tnum, N)
    d_ij = zeros(Tnum, N, N)
    b_ij = zeros(Tnum, N, N)
    q_ij = zeros(Tnum, Nbin)
    α_ij = zeros(Tnum, Nbin)
    ϵ_ij = zeros(Tnum, Nbin)

    c = map(Tnum, [0.16825491455291727, 10.296455569712531, 1.2476994961006087, 79.52567954035581, -1.5554400906899912, -100.34856564780112])
    zsum = sum(z)
    znorm = z ./ zsum

    pures = split_model(model)
    d_ii = @f(d)

    # Fill arrays
    for i in 1:N
        for j in i:N
            m_ii = pures[i].params.segment.values[1]
            m_jj = pures[j].params.segment.values[1]
            σ_ii = pures[i].params.sigma.values[1]
            σ_jj = pures[j].params.sigma.values[1]
            ϵ_ii[i] = pures[i].params.epsilon.values[1]
            if i == j
                m_ij[i,j] = m_ii
                σ_ij[i,j] = σ_ii
                d_ij[i,j] = d_ii[i]
            else
                avgm = (m_ii + m_jj) / Tnum(2)
                avgs = (σ_ii + σ_jj) / Tnum(2)
                di = d_ii[i] * znorm[i] + d_ii[j] * znorm[j]
                m_ij[i,j] = avgm
                m_ij[j,i] = avgm
                σ_ij[i,j] = avgs
                σ_ij[j,i] = avgs
                d_ij[i,j] = di
                d_ij[j,i] = di
            end
            bval = m_ij[i,j] * d_ij[i,j]^3
            b_ij[i,j] = bval
            b_ij[j,i] = bval
        end
    end

    b = sum(znorm[i] * b_ij[i,i] for i in 1:N)
    α_ii = [m_ij[i,i] * ϵ_ii[i] / T for i in 1:N]

    function Q_ii(α, b)
        (c[1]*log(b) + c[2])*α^2 + (c[3]*log(b) + c[4])*α + c[5]*log(b) + c[6]
    end

    # Calculate q_ij
    k = 1
    components = model.components
    for i in 1:N
        for j in i+1:N
            z_bin = [znorm[i]/(znorm[i]+znorm[j]), znorm[j]/(znorm[i]+znorm[j])]
            b_bin = z_bin[1]*b_ij[i,i] + z_bin[2]*b_ij[j,j]
            binary_model = Michelsen_GEPCSAFT([components[i], components[j]])
            gₑ = excess_gibbs_free_energy(binary_model.activity, V, T, z_bin) / (R̄ * T)
            b_sum = z_bin[1]*log(b_ij[i,i] / b_bin) + z_bin[2]*log(b_ij[j,j] / b_bin)
            α_sum = z_bin[1]*Q_ii(α_ii[i], b_ij[i,i]) + z_bin[2]*Q_ii(α_ii[j], b_ij[j,j])
            q_ij[k] = gₑ + b_sum + α_sum
            k += 1
        end
    end

    # Calculate α_ij and ϵ_ij
    k = 1
    for i in 1:N
        for j in i+1:N
            z_bin = [znorm[i]/(znorm[i]+znorm[j]), znorm[j]/(znorm[i]+znorm[j])]
            b_bin = z_bin[1]*b_ij[i,i] + z_bin[2]*b_ij[j,j]
            ϕ₁ = c[1]*log(b_bin) + c[2]
            ϕ₂ = c[3]*log(b_bin) + c[4]
            ϕ₃ = c[5]*log(b_bin) + c[6] - q_ij[k]
            Δ = ϕ₂^2 - 4*ϕ₁*ϕ₃
            Δ = real(Δ) # in case of negative (shouldn't happen)
            α_1 = (-ϕ₂ + sqrt(Δ)) / (2*ϕ₁)
            α_2 = (-ϕ₂ - sqrt(Δ)) / (2*ϕ₁)
            arithm = z_bin[1]*α_ii[i] + z_bin[2]*α_ii[j]
            if α_1 > zero(Tnum) && α_2 > zero(Tnum)
                α_ij[k] = abs(α_1 - arithm) < abs(α_2 - arithm) ? α_1 : α_2
            elseif α_1 > zero(Tnum)
                α_ij[k] = α_1
            else
                α_ij[k] = α_2
            end
            ϵ_ij[k] = T * α_ij[k] / m_ij[i,j]
            k += 1
        end
    end

    # Calculate m2ϵσ3₁, m2ϵσ3₂
    m2ϵσ3₁ = zero(Tnum)
    m2ϵσ3₂ = zero(Tnum)
    k = 1
    for i in 1:N
        for j in i+1:N
            m2ϵσ3₁ += znorm[i]*znorm[j]*m_ij[i,i]*m_ij[j,j]*ϵ_ij[k]*σ_ij[i,j]^3
            m2ϵσ3₂ += znorm[i]*znorm[j]*m_ij[i,j]*m_ij[j,j]*ϵ_ij[k]^2*σ_ij[i,j]^3
            k += 1
        end
    end
    return m2ϵσ3₁, m2ϵσ3₂
end

function I(model::Michelsen_GEPCSAFTModel, V, T, z, n , _data=@f(data))
    _,_,_,_,η,m̄ = _data
    if n == 1
        corr = PCSAFTconsts.corr1
    elseif n == 2
        corr = PCSAFTconsts.corr2
    end
    res = zero(η)
    @inbounds for i ∈ 1:7
        ii = i-1
        corr1,corr2,corr3 = corr[i]
        ki = corr1 + (m̄-1)/m̄*corr2 + (m̄-1)/m̄*(m̄-2)/m̄*corr3
        res +=ki*η^ii
    end
    return res
end

function Ī(model::Michelsen_GEPCSAFTModel, V, T, z, n , _data=@f(data))
    _,_,_,_,_,m̄ = _data
    if n == 1
        corr = PCSAFTconsts.corr1
    elseif n == 2
        corr = PCSAFTconsts.corr2
    end
    res = zero(m̄)
    @inbounds for i ∈ 1:7
        corr1,corr2,corr3 = corr[i]
        ki = corr1 + (m̄-1)/m̄*corr2 + (m̄-1)/m̄*(m̄-2)/m̄*corr3
        res +=ki
    end
    return res
end

function Ii(model::Michelsen_GEPCSAFTModel, V, T, z, n , _data=@f(data))
    m = model.params.segment.values
    if n == 1
        corr = PCSAFTconsts.corr1
    elseif n == 2
        corr = PCSAFTconsts.corr2
    end
    res = zero(m)
    @inbounds for i ∈ 1:7
        corr1,corr2,corr3 = corr[i]
        ki = @. corr1 + (m-1)/m*corr2 + (m-1)/m*(m-2)/m*corr3
        res +=ki
    end
    return res
end

function newtons_method(f, df, xn)
    return xn - f/df 
end