struct COFFEEParam <: EoSParam
    Mw::SingleParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    epsilon::PairParam{Float64}
    lambda_r::PairParam{Float64}
    lambda_a::PairParam{Float64}
    dipole::SingleParam{Float64}
    dipole2::SingleParam{Float64}
    shift::SingleParam{Float64}
end

abstract type COFFEEModel <: SAFTVRMieModel end

@newmodel COFFEE COFFEEModel COFFEEParam false
default_references(::Type{COFFEE}) = ["10.1063/1.5111364","10.1063/1.5136079"]
default_locations(::Type{COFFEE}) = ["COFFEE"]
default_ignore_missing_singleparams(::Type{COFFEE}) = ["dipole"]

"""
    COFFEEModel <: SAFTVRMieModel

    COFFEE(components;
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
## Description
Co-Oriented Fluid Functional Equation for Electrostatic interactions (COFFEE)
## References
1. Langenbach, K. (2017). Co-Oriented Fluid Functional Equation for Electrostatic interactions (COFFEE). Chemical Engineering Science, 174, 40-55 (2017). [doi:10.1016/j.ces.2017.08.025](https://doi.org/10.1016/j.ces.2017.08.025)
"""
COFFEE
export COFFEE

function transform_params(::Type{COFFEE},params,components)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    shift = params["shift"]
    shift.values .*= 1E-10
    sigma = sigma_LorentzBerthelot(sigma)
    epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], sigma)
    lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    params["sigma"] = sigma
    params["epsilon"] = epsilon
    params["lambda_a"] = lambda_a
    params["lambda_r"] = lambda_r
    μ = get!(params,"dipole") do
        SingleParam("dipole",components)
    end
    params["dipole2"] = SingleParam("Dipole squared",components, μ.^2 ./ k_B*1e-36*(1e-10*1e-3))
    return params
end

function a_res(model::COFFEEModel, V, T, z, _data = @f(data))
    return @f(a_hs,_data) + @f(a_dispchain,_data) + @f(a_ff,_data) + @f(a_nf,_data)
end

function a_ff(model ::COFFEEModel, V, T, z, _data=@f(data))
    dipole = pcp_dipole(model)
    if all(iszero,dipole)
        return zero(V+T+first(z))
    end
    a₂ = @f(a_2,_data)
    iszero(a₂) && return zero(a₂)
    a₃ = @f(a_3,_data)
    return a₂^2/(a₂-a₃)
end

function a_2(model ::COFFEEModel, V, T, z, _data=@f(data))
    (_,_,(_,_,_,η),_,_,_,_)= _data
    ∑z = sum(z)
    ρ = N_A*∑z/V
    _a_2 = zero(T+V+first(z))
    nc = length(model)
    ϵ = pcp_epsilon(model)
    σ = pcp_sigma(model)
    μ̄² = pcp_dipole2(model)
    @inbounds for i ∈ 1:nc
        _J2_ii = @f(J2,i,i,η,ϵ)
        zᵢ = z[i]
        μ̄²ᵢ = μ̄²[i]
        iszero(μ̄²ᵢ) && continue
        _a_2 +=zᵢ^2*μ̄²ᵢ^2/σ[i,i]^3*_J2_ii
        for j ∈ (i+1):nc
            μ̄²ⱼ = μ̄²[j]
            iszero(μ̄²ⱼ) && continue
            _J2_ij = @f(J2,i,j,η,ϵ)
            _a_2 += 2*zᵢ*z[j]*μ̄²ᵢ*μ̄²ⱼ/σ[i,j]^3*_J2_ij
        end
    end
    _a_2 *= -π*ρ/(T*T)/(∑z*∑z)
    return _a_2
end

function a_3(model ::COFFEEModel, V, T, z, _data=@f(data))
    ∑z = sum(z)
    ρ = N_A*∑z/V
    (_,_,(_,_,_,η),_,_,_,_)= _data
    _a_3 = zero(T+V+first(z))
    ϵ = pcp_epsilon(model)
    σ = pcp_sigma(model)
    μ̄² = pcp_dipole2(model)
    nc = length(model)

    @inbounds for i ∈ 1:nc
        zi,μ̄²i = z[i],μ̄²[i]
        iszero(μ̄²i) && continue
        _J3_iii = @f(J3,i,i,i,η)
        a_3_i = zi*μ̄²i/σ[i,i]
        _a_3 += a_3_i^3*_J3_iii
        for j ∈ i+1:nc
            zj,μ̄²j = z[j],μ̄²[j]
            iszero(μ̄²j) && continue
            σij⁻¹ = 1/σ[i,j]
            a_3_iij = zi*μ̄²i*σij⁻¹
            a_3_ijj = zj*μ̄²j*σij⁻¹
            a_3_j = zj*μ̄²j/σ[j,j]
            _J3_iij = @f(J3,i,i,j,η)
            _J3_ijj = @f(J3,i,j,j,η)
            _a_3 += 3*a_3_iij*a_3_ijj*(a_3_i*_J3_iij + a_3_j*_J3_ijj)
            for k ∈ j+1:nc
                zk,μ̄²k = z[k],μ̄²[k]
                iszero(μ̄²k) && continue
                _J3_ijk = @f(J3,i,j,k,η)
                _a_3 += 6*zi*zj*zk*μ̄²i*μ̄²j*μ̄²k*σij⁻¹/σ[i,k]/σ[j,k]*_J3_ijk
            end
        end
    end
    _a_3 *= -4π^2/3*ρ^2/(T*T*T)/(∑z*∑z*∑z)
    return _a_3
end

function J2(model::COFFEEModel, V, T, z, i, j, η = @f(ζ,3),ϵ = pcp_epsilon(model))
    b2 = COFFEEconsts.corr_b2
    c2 = COFFEEconsts.corr_c2
    ϵT = ϵ[i,j]/T
    c = b2 .+ c2 .* ϵT
    return evalpoly(η,c)
end

function J3(model::COFFEEModel, V, T, z, i, j, k, η = @f(ζ,3),ϵ = pcp_epsilon(model))
    b3 = COFFEEconsts.corr_b3
    c3 = COFFEEconsts.corr_c3
    ϵT = ϵ[i,j]/T
    c = b3 .+ c3 .* ϵT
    return evalpoly(η,c)
end

function a_nf(model ::COFFEEModel, V, T, z, _data=@f(data))
    (_,_,(_,_,_,η),_,_,_,_)= _data

    ϵ = model.params.epsilon[1]
    σ = diagvalues(model.params.sigma.values)
    d = model.params.shift[1] / σ[1]
    μ² = model.params.dipole2[1]./ϵ/σ[1]^3
    μ = sqrt(μ²)

    ρ̄ = N_A/V*sum(z.*σ.^3)
    g_hs = (1-η/2)/(1-η)^3
    T̄ = T/ϵ

    _Iμμ = Iμμ(ρ̄,T̄,μ)

    Q = ∫∫∫Odξ₁dξ₂dγ12(ρ̄,T̄,μ²,d,_Iμμ)

    return 19*π/12*ρ̄*g_hs*log(4π/Q)
end

function Iμμ(ρ̄,T̄,μ)
    a = COFFEEconsts.corr_I

    return a[1] + a[2]/T̄ + a[3]/T̄^2 + ρ̄ * (
           a[4] + a[5]*μ + a[6]/T̄ + a[7]/T̄^2 + a[8]*ρ̄ + a[9]*ρ̄^2 + a[10]*ρ̄*μ/T̄
    )
end

function ∫odr(d,ξ1,ξ2,γ12)
    if d==0.0
        cosΘ = 2*ξ1*ξ2 - √(1-ξ1^2)*√(1-ξ2^2)*cos(γ12)
        _I = -cosΘ*5/18
    else
        A12 = ξ1*ξ2 + √(1-ξ1^2)*√(1-ξ2^2)*cos(γ12)
        a = 3*ξ1*ξ2-A12
        b = d*(ξ1-ξ2)*(A12-3)
        c = -d^2*(A12^2-4A12+3)
        f = 2d*(ξ1-ξ2)
        f2 = f*f
        f3 = f2*f

        g = -2d^2*(A12-1)
        g2 = g*g

        _I = -(1/(3*(f2 - 4*g)^2))*((1/((1 + f + g)^(
            3/2)))*(2*c*(2 + f)*(-8 - 8*f + f2 - 12*g) +
            2*b*(3*f3 + 8*g2 + 2*f2*(6 + g) + 4*f*(2 + 3*g)) -
            2*a*(3*f3 + 8*g + 4*f*g*(3 + 2*g) + 2*f2*(1 + 6*g))) + (
            1/((9 + 6*f + 4*g)^(3/2))) *
            4*(-4*c*(3 + f)*(-12*f + f2 - 6*(3 + 2*g)) -
            2*b*(9*f3 + 16*g2 + 18*f*(3 + 2*g) + f2*(54 + 4*g)) +
            a*(27*f3 + 108*g + 9*f2*(3 + 8*g) + 4*f*g*(27 + 8*g))))
    end
    return _I
end

function ∫∫∫Odξ₁dξ₂dγ12(ρ̄,T̄,μ²,d, _Iμμ = Iμμ(ρ̄,T̄,μ))
    I(x) = exp(-24/19*μ²/T̄*_Iμμ*∫odr(d,x[1],x[2],x[3]))

    _I = ∫∫∫dξ₁dξ₂dγ12(I)

    return _I
end

function ∫∫∫Oodξ₁dξ₂dγ12(ρ̄,T̄,μ²,d,_Iμμ,Q)
    I(x) = begin
        _I = ∫odr(d,x[1],x[2],x[3])
        return exp(-24/19*μ²/T̄*_Iμμ*_I)*_I
    end

    _I = ∫∫∫dξ₁dξ₂dγ12(I)

    return _I/Q
end


function ∫∫∫OlnOdξ₁dξ₂dγ12(ρ̄,T̄,μ²,d,_Iμμ,Q)
    I(x) = begin
        _I = -24/19*μ²/T̄*_Iμμ*∫odr(d,x[1],x[2],x[3])
        return _I*exp(_I)
    end

    _I = ∫∫∫dξ₁dξ₂dγ12(I)

    return _I/Q
end

function ∫∫∫dξ₁dξ₂dγ12(I::Function)
    ξ1,w1 = COFFEEconsts.xi_quadrature
    ξ2,w2 = COFFEEconsts.xi_quadrature
    γ12,w3 = COFFEEconsts.gamma_quadrature

    _I = 0.

    nξ = length(ξ1)
    nγ = length(γ12)

    for i in 1:nξ
        for j in 1:nξ
            for k in 1:nγ
                _I += w1[i]*w2[j]*w3[k]*I([ξ1[i],ξ2[j],γ12[k]])
            end
        end
    end
    return _I
end

#=
function Xij(model::COFFEEModel, V, T, z, _data=@f(data))
    nc = length(model)
    Χᵢⱼ = similar(Base.promote_eltype(model,V,T,z),nc,nc)
    ξᵢⱼ = similar(Χᵢⱼ)
    K = similar(Χᵢⱼ)
    ∑z = sum(z)
    σ = pcp_sigma(model)
    Χᵢⱼ .= 0
    #filling initial values
    for i in 1:nc
        Χᵢⱼ[i,i] = z[i]*z[i]/(∑z*∑z)
        ξᵢⱼ[i,i] = COFFEE_ξᵢⱼ(z,σ,i,i)
        for j in i:nc
            Χᵢⱼ[i,j] = z[i]*z[j]/(∑z*∑z)
            ξᵢⱼ[i,j] = COFFEE_ξᵢⱼ(z,σ,i,j)
            K[i,j] = (19/6)*σ[i,j]^3 * 1 #TODO: put gij *  ∫Oijdωidωj here
        end
    end
    Χnew = copy(Χᵢⱼ)
    tol = one(eltype(Χᵢⱼ))
    while tol > 1e-10
        Xij_fixpoint(model::COFFEEModel,V,T,z,data,Χnew,Χ,ξ,K)
        tol = Solvers.dnorm(Χᵢⱼ,Χnew)
        Χᵢⱼ .= Χnew
    end
    return Χᵢⱼ
end

function COFFEE_Δij(model::COFFEEModel,V,T,z,data,Χ,i,j)
    return 1.0 #TODO: fill
end

function Xij_fixpoint(model::COFFEEModel,V,T,z,data,Χnew,Χ,ξ,K)
    ∑z = sum(z)
    for i in 1:nc
        #=
        local mass balance
        =#
        Χii = zero(eltype(Χnew))
        σii = σ[i,i]
        zᵢ = z[i]
        for j in 1:nc
            if i !== j
                σjj = σ[j,j]
                Χii += z[i]*σjj*σjj*σjj*(1 - Χ[i,j])
            end
        end
        Χii = Χii/(zᵢ*σii*σii*σii)

        if Χii < 0
            Χii = Χ[i,i] / 2
        elseif Χii > 1
            Χii =  1 - 0.5*Χ[i,i]
        end

        Χnew[i,i] = Χii

        #=
        Energy minimization constraint
        =#
        for j in i+1:nc
            Δᵢⱼ = @f(COFFEE_Δij,data,Χ,i,j)
            Χᵢⱼ = Χ[i,j]
            ξᵢⱼ = ξ[i,j]
            zⱼ = z[j]
            Χij = ∑z/(zᵢ*(1 + 1/ξᵢⱼ) + zⱼ*(1 + ξᵢⱼ)) * K[i,j] * log(ξᵢⱼ + Χᵢⱼ*(zᵢ - zⱼ*ξᵢⱼ)/∑z) / Δᵢⱼ

            if Χij < 0
                Χij = Χᵢⱼ / 2
            elseif Χij > 1
                Χij =  1 - 0.5*Χᵢⱼ
            end

            Χnew[i,j] = Χij

        end
    end
    return Χnew
end

function COFFEE_ξᵢⱼ(z,σ,i,j)
    k = length(z)
    ∑1 = zero(Base.promote_eltype(z,σ))
    ∑2 = zero(Base.promote_eltype(z,σ))
    for k in 1:length(z)
        σjk = σ[j,k]
        σik = σ[i,k]
        if k != i
            ∑1 += z[k]*σjk^3
        end
        if k != j
            ∑2 += z[k]*σik^3
        end
    end
    return ∑1/∑2
end
=#
const COFFEEconsts = (
    corr_I = (1.4713,-0.7842,0.6010,1.2336,-0.6978,-1.6837,0.4262,-0.4824,0.3502,0.2707),

    corr_b2 =
    ( 0.106652,
    -0.0710423,
      1.216337,
     0.1544786,
     -3.147176),

    corr_b3 =
    (-6.267389e-3,
     -12.49881e-3,
      40.96034e-3,
    -0.6194877e-3),

    corr_c2 =
    (0.3322687,
     -1.511604,
      1.201297,
      0.0,
      0.0),

    corr_c3 =
    (-13.53321e-3,
       37.2739e-3,
    -0.1944546e-3,
     0.0),

    xi_quadrature = ((
        -0.9955569697904981, -0.9766639214595175, -0.9429745712289743, -0.8949919978782753, -0.833442628760834, -0.7592592630373577, -0.6735663684734684, -0.577662930241223, -0.473002731445715, -0.36117230580938786, -0.24386688372098844, -0.1228646926107104, 0.0, 0.1228646926107104, 0.24386688372098844, 0.36117230580938786, 0.473002731445715, 0.577662930241223, 0.6735663684734684, 0.7592592630373577, 0.833442628760834, 0.8949919978782753, 0.9429745712289743, 0.9766639214595175, 0.9955569697904981),
                     (0.01139379850102625, 0.026354986615031935, 0.04093915670130625, 0.054904695975835305, 0.06803833381235695, 0.08014070033500113, 0.09102826198296363, 0.10053594906705073, 0.10851962447426364, 0.11485825914571161, 0.11945576353578476, 0.12224244299031012, 0.12317605372671545, 0.12224244299031012, 0.11945576353578476, 0.11485825914571161, 0.10851962447426364, 0.10053594906705073, 0.09102826198296363, 0.08014070033500113, 0.06803833381235695, 0.054904695975835305, 0.04093915670130625, 0.026354986615031935, 0.01139379850102625)),
    gamma_quadrature = ((0.010793571159988413, 0.0565927642933524, 0.13786183772187116, 0.2527144697052986, 0.3984609095574594, 0.5716855414653129, 0.768328316649816, 0.9837801753398278, 1.2129911485257514, 1.450588748496007, 1.6910039050937862, 1.9286015050640417, 2.1578124782499652, 2.373264336939977, 2.5699071121244805, 2.743131744032334, 2.888878183884495, 3.003730815867922, 3.0849998892964408, 3.130799082429805),
                        (0.027668017714319135, 0.06377657679306871, 0.09844502331593069, 0.13081079977613538, 0.16011145779868483, 0.18565953665239512, 0.20685602955658772, 0.22320404656916054, 0.2343203792081908, 0.239944459410423, 0.239944459410423, 0.2343203792081908, 0.22320404656916054, 0.20685602955658772, 0.18565953665239512, 0.16011145779868483, 0.13081079977613538, 0.09844502331593069, 0.06377657679306871, 0.027668017714319135)),
)