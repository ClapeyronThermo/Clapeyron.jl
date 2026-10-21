struct SAFTVRQMieParam <: EoSParam
    Mw::PairParam{Float64}
    segment::SingleParam{Float64}
    sigma::PairParam{Float64}
    lambda_a::PairParam{Float64}
    lambda_r::PairParam{Float64}
    epsilon::PairParam{Float64}
end

abstract type SAFTVRQMieModel <: SAFTVRMieModel end
#SAFT-VRQ-Mie does not have sites.

struct SAFTVRQMie{I} <: SAFTVRQMieModel
    components::Vector{String}
    params::SAFTVRQMieParam
    idealmodel::I
    fh_order::Symbol
    references::Array{String,1}
end

function SAFTVRQMie(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    fh_order = :fh2,
    verbose = false)

    fh_order = Symbol(fh_order)

    MODEL,PARAM = SAFTVRQMie,SAFTVRQMieParam
    locations = default_locations_vrq(fh_order)
    _components = format_components(components)
    params_in = getparams(_components, locations,return_sites = false)
    params_out = transform_params(MODEL,params_in)
    pkgparam = build_eosparam(PARAM,params_out)
    init_idealmodel = init_model(idealmodel,_components,ideal_userlocations,verbose)
    model = SAFTVRQMie(_components,pkgparam,init_idealmodel,fh_order,default_references(MODEL))
    set_reference_state!(model,reference_state;verbose)
    return model
end

default_references(::Type{SAFTVRQMie}) = ["10.1063/1.5111364","10.1063/1.5136079"]

default_locations(::Type{SAFTVRQMie}) = default_locations_vrq(:fh2)
default_locations(model::SAFTVRQMieModel) = default_locations_vrq(model.fh_order)

function default_locations_vrq(fh_order)
    if fh_order == :fh2
        locations = ["SAFT/SAFTVRQMie/fh2"]
    elseif fh_order == :fh1
        locations = ["SAFT/SAFTVRQMie/fh1"]
    else
        error("Unknown Feynman-Hibbs order: $fh_order")
    end
end

function transform_params(::Type{SAFTVRQMie},params)
    Mw = params["Mw"]
    Mw .*= 1E-3
    mw_mix(mi,mj,k) = mix_powmean(mi,mj,k,-1) #mij = 0.5/(mi^-1 + mj^-1)
    Mw = kij_mix(mw_mix,Mw)
    k = get(params,"k",nothing)
    l = get(params,"l",nothing)
    params["sigma"].values .*= 1E-10
    sigma = sigma_LorentzBerthelot(params["sigma"],l)
    epsilon = epsilon_HudsenMcCoubrey(params["epsilon"], sigma)
    lambda_a = lambda_LorentzBerthelot(params["lambda_a"])
    lambda_r = lambda_LorentzBerthelot(params["lambda_r"])
    if k !== nothing
        epsilon .= epsilon .* (1 .- k)
    end
    params["Mw"] = Mw
    params["sigma"] = sigma
    params["epsilon"] = epsilon
    params["lambda_a"] = lambda_a
    params["lambda_r"] = lambda_r
    return params
end

function show_info(io,model::SAFTVRQMie)
    fh = model.fh_order
    println(io)
    if fh == :fh1
        print(io,"Feynman-Hibbs Perturbation: 1st order")
    elseif fh == :fh2
        print(io,"Feynman-Hibbs Perturbation: 2nd order")
    end
end

"""
    SAFTVRQMieModel <: SAFTVRMieModel

    SAFTVRQMie(components;
    idealmodel = BasicIdeal,
    userlocations = String[],
    ideal_userlocations = String[],
    reference_state = nothing,
    fh_order = :fh2,
    verbose = false)

## Input parameters
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Single Parameter (`Float64`) - Segment Diameter `[Å]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon`: Single Parameter (`Float64`) - Reduced dispersion energy `[K]`

## Model Parameters
- `Mw`: Pair Parameter (`Float64`) - Mixed Molecular Weight `[g·mol⁻¹]`
- `segment`: Single Parameter (`Float64`) - Number of segments (no units)
- `sigma`: Pair Parameter (`Float64`) - Mixed segment Diameter `[m]`
- `lambda_a`: Pair Parameter (`Float64`) - Atractive range parameter (no units)
- `lambda_r`: Pair Parameter (`Float64`) - Repulsive range parameter (no units)
- `epsilon`: Pair Parameter (`Float64`) - Mixed reduced dispersion energy`[K]`

## Input models
- `idealmodel`: Ideal Model

## Description

Quantum-Corrected SAFT-VR Mie. In particular,The Feynman–Hibbs correction order can be modified by passing the `fh_order` keyword argument. The default is 2nd order (`:fh2`), but 1st order (`:fh1`) is also available.

## References
1. Aasen, A., Hammer, M., Ervik, Å., Müller, E. A., & Wilhelmsen, Ø. (2019). Equation of state and force fields for Feynman–Hibbs-corrected Mie fluids. I. Application to pure helium, neon, hydrogen, and deuterium. The Journal of Chemical Physics, 151(6), 064508. [doi:10.1063/1.5111364](https://doi.org/10.1063/1.5111364)
2. Aasen, A., Hammer, M., Müller, E. A., & Wilhelmsen, Ø. (2020). Equation of state and force fields for Feynman-Hibbs-corrected Mie fluids. II. Application to mixtures of helium, neon, hydrogen, and deuterium. The Journal of Chemical Physics, 152(7), 074507. [doi:10.1063/1.5136079](https://doi.org/10.1063/1.5136079)
"""
SAFTVRQMie

export SAFTVRQMie

mw(model::SAFTVRQMieModel) = diagvalues(model.params.Mw) .* 1e3

function a_mono(model::SAFTVRQMieModel, V, T, z,_data = @f(data))
    ahs = @f(a_hs_eff,_data)
    adisp = @f(a_disp,_data)
    return ahs + adisp
end

function a_res(model::SAFTVRQMieModel, V, T, z)
    @f(a_mono)
end

function data(model::SAFTVRQMieModel, V, T, z)
    m̄ = dot(z,model.params.segment.values)
    _σeff = @f(σeff)
    _ϵff = @f(ϵeff)
    _d = @f(d,_σeff) #d here is a (comp x comp) matrix, instead of a (comp) vector, as all safts
    ζi = @f(ζ0123,diagvalues(_d))
    _ζst = @f(ζst,_σeff)
    ζₓ = @f(ζ_X,_d)
    _ρ_S = @f(ρ_S,m̄)
    σ3x = _ζst/(_ρ_S*π/6)
    vrdata = (_d,_ρ_S,ζi,ζₓ,_ζst,σ3x,m̄)
    return (_σeff,_ϵff,vrdata)
end


function ζ_X(model::SAFTVRQMieModel, V, T, z,_d = @f(d))
    m = model.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    ρS = N_A/V*m̄
    comps = 1:length(z)
    ∑xixjdij³ = zero(first(z))
    for i ∈ comps
        x_Si = z[i]*m[i]*m̄inv
        di =_d[i,i]
        r1 = x_Si*x_Si*(di)^3
        ∑xixjdij³ += r1
        for j ∈ 1:(i-1)
            x_Sj = z[j]*m[j]*m̄inv
            dij = _d[i,j]
            r1 = x_Si*x_Sj*dij^3
            ∑xixjdij³ += 2*r1
        end
    end
    ζₓ = ∑xixjdij³*ρS* π/6
    return ζₓ
end

vrq_Q1(λ) = λ*(λ-1)
vrq_Q2(λ,second_order) = second_order ? 1/2*(λ+2)*(λ+1)*λ*(λ-1) : zero(λ/2)
Q1(model::SAFTVRQMieModel, V, T, z, λ) = vrq_Q1(λ)
Q2(model::SAFTVRQMieModel, V, T, z, λ) = vrq_Q2(λ,model.fh_order === :fh2)

function vrq_u(x, parameters)
    Q1λr, Q1λa, Q2λr, Q2λa, D, λr, λa, C, _, second_order = parameters
    u_over_C = (x^-λr - x^-λa) +
                D*(Q1λr*x^-(λr+2) - Q1λa*x^-(λa+2))
    if second_order
        u_over_C += D*D*(Q2λr*x^-(λr+4) - Q2λa*x^-(λa+4))
    end
    return u_over_C * C
end

function deff(σ, x_eff, params)
    function deff_obj(x, parameters)
        β = parameters[end-1]
        return exp(-β * vrq_u(x,parameters))
    end
    fgh(xx,p) = Solvers.to_halley(Base.Fix2(deff_obj,p),xx)
    prob  = Roots.ZeroProblem(fgh, one(1.0*primalval(x_eff)))
    x_min = roots_solve_ad(prob, Roots.Halley(), params)
    f̄     = Base.Fix2(deff_obj, params)
    return σ * (x_eff - Solvers.integral21(f̄, x_min, x_eff))
end

#compatibility
function deff(model::SAFTVRQMieModel, V, T, z, _data)
    ϵ,σ,λr,λa,Mw,σeff = _data
    Q1λr,Q1λa,Q2λr,Q2λa = @f(Q1,λr),@f(Q1,λa),@f(Q2,λr),@f(Q2,λa)
    D = ħ^2/(12*k_B*T*Mw/N_A*σ^2)
    C = @f(Cλ,λa,λr)
    β = ϵ/T
    is_second_order_fh = model.fh_order == :fh2
    params = (Q1λr,Q1λa,Q2λr,Q2λa,D,λr,λa,C,β,is_second_order_fh)
    return deff(σ, σeff/σ, params)
end

function d(model::SAFTVRQMieModel, V, T, z,_σeff = @f(σeff))
    _ϵ = model.params.epsilon.values
    _σ = model.params.sigma.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    _Mw = model.params.Mw.values
    n = length(model)
    _d = zeros(Base.promote_eltype(model,T,1.0),n,n)
    is_second_order_fh = model.fh_order == :fh2
    for i ∈ 1:n
        for j ∈ 1:i
            σeffij = _σeff[i,j]
            λaij,λrij,σij,Mwij = _λa[i,j],_λr[i,j],_σ[i,j],_Mw[i,j]
            Q1λrij,Q1λaij,Q2λrij,Q2λaij = @f(Q1,λrij),@f(Q1,λaij),@f(Q2,λrij),@f(Q2,λaij)
            Dij = ħ^2/(12*k_B*T*Mwij/N_A*σij^2)
            Cij = @f(Cλ,λaij,λrij)
            βij = _ϵ[i,j]/T
            params_ij = (Q1λrij,Q1λaij,Q2λrij,Q2λaij,Dij,λrij,λaij,Cij,βij,is_second_order_fh)  
            dij = deff(σij, σeffij/σij, params_ij)
            _d[i,j] = dij
            if i != j
                _d[j,i] = dij
            end
        end
    end
    _d
end

#specialization, and solves ambiguity
function d(model::SAFTVRQMieModel, V, T, z::SingleComp,_σeff = @f(σeff))
    _ϵ = model.params.epsilon[1]
    _σ = model.params.sigma[1]
    _λr = model.params.lambda_r[1]
    _λa = model.params.lambda_a[1]
    _Mwij = model.params.Mw[1]
    _data = (_ϵ,_σ,_λr,_λa,_Mwij,_σeff[1,1])
    d_ij = @f(deff,_data)
    return SA[d_ij;;]
end

function vrq_σeff_obj(x, parameters)
    Q1λr, Q1λa, Q2λr, Q2λa, D, λr, λa, second_order = parameters
    xΔλ = x^(λr-λa)
    x2  = x*x
    x4  = x2*x2
    f   = D*Q1λr*x2 + x4 - xΔλ*(D*Q1λa*x2 + x4)
    if second_order
        f += D*D*(Q2λr - xΔλ*Q2λa)
    end
    return f
end

function σeff(model::SAFTVRQMieModel, V, T, z)
    _σ = model.params.sigma.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    _Mw = model.params.Mw.values
    #Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)
    #f(x) = Dij^2*@f(Q2,λr[i,j])+Dij*@f(Q1,λr[i,j])*x^2+x^4-x^(λr[i,j]-λa[i,j])*(Dij^2*@f(Q2,λa[i,j])+Dij*@f(Q1,λa[i,j])*x^2+x^4)
    #g(x) = 2*Dij*@f(Q1,λr[i,j])*x+4x^3-x^(λr[i,j]-λa[i,j]-1)*((λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+(λr[i,j]-λa[i,j]+2)*Dij*@f(Q1,λa[i,j])*x^2+(λr[i,j]-λa[i,j]+4)*x^4)
    #h(x) = 2*Dij*@f(Q1,λr[i,j])+12x^2-x^(λr[i,j]-λa[i,j]-2)*((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+2)+2*(λr[i,j]-λa[i,j]+2))*Dij*@f(Q1,λa[i,j])*x^2+((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+4)+4*(λr[i,j]-λa[i,j]+4))*x^4)


    fgh(x,p) = Solvers.to_halley(Base.Fix2(vrq_σeff_obj,p),x)
    prob = Roots.ZeroProblem(fgh,one(1.0*primalval(T)))
    is_second_order_fh = model.fh_order == :fh2

    _σeff = zeros(Base.promote_eltype(model,T,1.0),size(_σ))
    nc = length(model)

    #=
    roots_solve_ad performs a Roots.jl 1-D rootfinding, followed by solving the derivative of the solution respect to the parameters, via the implicit function theorem.
    =#
    for i in 1:nc
        for j in 1:i
            λaij,λrij,σij,Mwij = _λa[i,j],_λr[i,j],_σ[i,j],_Mw[i,j]
            Q1λrij,Q1λaij,Q2λrij,Q2λaij = @f(Q1,λrij),@f(Q1,λaij),@f(Q2,λrij),@f(Q2,λaij)
            Dij = ħ^2/(12*k_B*T*Mwij/N_A*σij^2)
            params_ij = (Q1λrij,Q1λaij,Q2λrij,Q2λaij,Dij,λrij,λaij,is_second_order_fh)
            σeff_ij = σij*roots_solve_ad(prob,Roots.Halley(),params_ij)
            _σeff[i,j] = σeff_ij
            if i != j
                _σeff[j,i] = σeff_ij
            end
        end
    end
    return _σeff
end

function σeff(model::SAFTVRQMieModel, V, T, z::SingleComp)
    _σ  = model.params.sigma[1]
    _λr = model.params.lambda_r[1]
    _λa = model.params.lambda_a[1]
    _Mw = model.params.Mw[1]
    D   = ħ^2/(12*k_B*T*_Mw/N_A*_σ^2)
    is_second_order_fh = model.fh_order == :fh2
    fgh(x,p) = Solvers.to_halley(Base.Fix2(vrq_σeff_obj,p), x)
    prob = Roots.ZeroProblem(fgh, one(1.0*primalval(T)))
    params = (@f(Q1,_λr), @f(Q1,_λa), @f(Q2,_λr), @f(Q2,_λa),
              D, _λr, _λa, is_second_order_fh)

    σeff_ij = _σ * roots_solve_ad(prob, Roots.Halley(), params)
    return SA[σeff_ij;;]
end

function vrq_ϵeff_obj(x, parameters)
    Q1λr, Q1λa, Q2λr, Q2λa, D, λr, λa, _, _ ,second_order = parameters
    x2  = x*x
    x4  = x2*x2
    xΔλ = x^(λr-λa)
    f   = (λr+2)*D*Q1λr*x2 + λr*x4 -
            xΔλ*((λa+2)*D*Q1λa*x2 + λa*x4)
    if second_order
        f += D*D*((λr+4)*Q2λr - xΔλ*(λa+4)*Q2λa)
    end
    return f
end

function ϵeff(model::SAFTVRQMieModel, V, T, z)
    _ϵ = model.params.epsilon.values
    _σ = model.params.sigma.values
    _λr = model.params.lambda_r.values
    _λa = model.params.lambda_a.values
    _Mw = model.params.Mw.values
    #Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)
    #f(x) = (λr[i,j]+4)*Dij^2*@f(Q2,λr[i,j])+(λr[i,j]+2)*Dij*@f(Q1,λr[i,j])*x^2+λr[i,j]*x^4-x^(λr[i,j]-λa[i,j])*((λa[i,j]+4)*Dij^2*@f(Q2,λa[i,j])+(λa[i,j]+2)*Dij*@f(Q1,λa[i,j])*x^2+λa[i,j]*x^4)
    #g(x) = 2*(λr[i,j]+2)*Dij*@f(Q1,λr[i,j])*x+4λr[i,j]*x^3-x^(λr[i,j]-λa[i,j]-1)*((λa[i,j]+4)*(λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+(λa[i,j]+2)*(λr[i,j]-λa[i,j]+2)*Dij*@f(Q1,λa[i,j])*x^2+λa[i,j]*(λr[i,j]-λa[i,j]+4)*x^4)
    #h(x) = 2*(λr[i,j]+2)*Dij*@f(Q1,λr[i,j])+12λr[i,j]*x^2-x^(λr[i,j]-λa[i,j]-2)*((λa[i,j]+4)*(λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j])*Dij^2*@f(Q2,λa[i,j])+(λa[i,j]+2)*((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+2)+2*(λr[i,j]-λa[i,j]+2))*Dij*@f(Q1,λa[i,j])*x^2+λa[i,j]*((λr[i,j]-λa[i,j]-1)*(λr[i,j]-λa[i,j]+4)+4*(λr[i,j]-λa[i,j]+4))*x^4)
    #σ_min = @f(Halley,f,g,h,(λr[i,j]/λa[i,j])^(1/(λr[i,j]-λa[i,j])))
    #return -ϵ[i,j]*u(σ_min)

    is_second_order_fh = model.fh_order == :fh2
    _ϵeff = zeros(Base.promote_eltype(model,T),size(_σ))
    nc = length(model)
    fgh(x,p) = Solvers.to_halley(Base.Fix2(vrq_ϵeff_obj,p),x)
    #=
    roots_solve_ad performs a Roots.jl 1-D rootfinding, followed by solving the derivative of the solution respect to the parameters, via the implicit function theorem.
    =#
    for i in 1:nc
        for j in 1:i
            λaij,λrij,σij,Mwij = _λa[i,j],_λr[i,j],_σ[i,j],_Mw[i,j]
            Q1λrij,Q1λaij,Q2λrij,Q2λaij = @f(Q1,λrij),@f(Q1,λaij),@f(Q2,λrij),@f(Q2,λaij)
            Dij = ħ^2/(12*k_B*T*Mwij/N_A*σij^2)
            Cij = @f(Cλ,λaij,λrij)
            β = Cij #not used, but defined to make the tuple complete.
            params_ij = (Q1λrij,Q1λaij,Q2λrij,Q2λaij,Dij,λrij,λaij,Cij,β,is_second_order_fh)
            x0_ij = (λrij/λaij)^(1/(λrij-λaij))
            prob_ij = Roots.ZeroProblem(fgh,x0_ij)
            xmin_ij = roots_solve_ad(prob_ij,Roots.Halley(),params_ij)
            uij = vrq_u(xmin_ij,params_ij)
            ϵeff_ij = -_ϵ[i,j]*uij
            _ϵeff[i,j] = ϵeff_ij
            if i != j
                _ϵeff[j,i] = ϵeff_ij
            end
        end
    end
    return _ϵeff
end

function ϵeff(model::SAFTVRQMieModel, V, T, z::SingleComp)
    _ε  = model.params.epsilon[1]
    _σ  = model.params.sigma[1]
    _λr = model.params.lambda_r[1]
    _λa = model.params.lambda_a[1]
    _Mw = model.params.Mw[1]
    D   = ħ^2/(12*k_B*T*_Mw/N_A*_σ^2)
    is_second_order_fh = model.fh_order == :fh2
    C = @f(Cλ,_λa, _λr)
    params = (@f(Q1,_λr), @f(Q1,_λa), @f(Q2,_λr), @f(Q2,_λa),
              D, _λr, _λa, C, C, is_second_order_fh)

    fgh(x,p) = Solvers.to_halley(Base.Fix2(vrq_ϵeff_obj,p), x)
    x0_ij = (_λr/_λa)^(1/(_λr-_λa))
    prob  = Roots.ZeroProblem(fgh, x0_ij)
    xmin  = roots_solve_ad(prob, Roots.Halley(), params)
    ϵeff_ij = -_ε * vrq_u(xmin, params)
    return SA[ϵeff_ij;;]
end

#=
Hard-sphere reference of SAFT-VRQ Mie for mixtures: non-additive hard spheres, Santos, López de Haro & Yuste
(J. Chem. Phys. 122, 024514 (2005)), as used in Aasen et al., J. Chem. Phys. 152, 074507 (2020), eqs. 14-24:

    d_na³ = Σ x_i d_ii³,   η_na = π ρ d_na³ / 6
    B̄2 = 4 Σ_ij x_i x_j d_ij³                          (exact B2 of the non-additive mixture)
    B̄3 = Σ_ijk x_i x_j x_k B̄_ijk,  B̄_ijk = 4/3 (c_k;ij d_ij³ + c_j;ik d_ik³ + c_i;jk d_jk³)
    c_k;ij = δ_k;ij³ + 3/2 δ_k;ij² δ_i;jk δ_j;ik / d_ij,   δ_k;ij = max(d_ik + d_jk - d_ij, 0)
    A1 = (10 d_na³ B̄2 - 4 B̄3)/(6 d_na⁶),   A2 = (B̄3 - d_na³ B̄2)/(6 d_na⁶)
    a_HS = -ln(1-η_na) A1 + a_CS(η_na) A2                (per particle/segment)

Sums use segment fractions x_S,i = z_i m_i / Σ z_j m_j (= mole fractions when m_i = 1, as in the paper);
the result is per mole of mixture (× m̄/Σz), consistent with a_disp.
=#
function a_hs_eff(model::SAFTVRQMieModel, V, T, z, _data = @f(data))
    _σeff, _ϵff, vrdata = _data
    _d, _ρ_S, _, _, _, _, m̄ = vrdata          # m̄ = Σ z_i m_i (not normalised)
    m    = model.params.segment.values
    nc   = length(z)
    Σz   = sum(z)
    _0   = zero(V + T + first(z))

    d3_na = _0
    for i in 1:nc
        xSi = z[i]*m[i]
        d3_na += xSi*_d[i,i]^3
    end
    d3_na /= m̄
    η_na = _ρ_S*(π/6)*d3_na
    #N_A/V*m̄
    B̄₂ = _0
    for i in 1:nc
        xSi = z[i]*m[i]
        B̄₂ += xSi*xSi*_d[i,i]^3
        for j in 1:(i-1)
            xSj = z[j]*m[j]
            B̄₂ += 2*xSi*xSj*_d[i,j]^3
        end
    end    
    B̄₂ *= 4/(m̄*m̄)

    B̄₃ = _0
    for i in 1:nc, j in i:nc, k in j:nc
        xSi,xSj,xSk = z[i]*m[i],z[j]*m[j],z[k]*m[k]
        dij, dik, djk = _d[i,j], _d[i,k], _d[j,k]
        δk = max(dik + djk - dij, _0)      # δ_k;ij
        δi = max(dij + dik - djk, _0)      # δ_i;jk
        δj = max(dij + djk - dik, _0)      # δ_j;ik
        c_kij = δk^3 + 1.5*δk*δk*δi*δj/dij
        c_ijk = δi^3 + 1.5*δi*δi*δj*δk/djk
        c_jik = δj^3 + 1.5*δj*δj*δi*δk/dik
        mult = (i == j == k) ? 1 : ((i == j || j == k) ? 3 : 6)
        B̄₃ += mult*xSi*xSj*xSk*(c_kij*dij^3 + c_ijk*djk^3 + c_jik*dik^3)
    end
    B̄₃ *= 4/(3*m̄*m̄*m̄)

    A1 = (10*d3_na*B̄₂ - 4*B̄₃)/(6*d3_na^2)
    A2 = (B̄₃ - d3_na*B̄₂)/(6*d3_na^2)
    a_CS = η_na*(4 - 3*η_na)/(1 - η_na)^2
    a_seg = -log1p(-η_na)*A1 + a_CS*A2

    return a_seg*m̄/Σz
end
#=
function a_1(model::SAFTVRQMieModel, V, T, z, i, j,_data = @f(data))
    _σeff,_ϵff,_d,_ρ_S,ζi,ζₓ,_ζst = _data
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)
    x_0ij = @f(x_0,i,j)
    dij = 0.5*(_d[i] + _d[j])
    x_0effij = _σeff[i,j]/dij
    return 2*π*ϵ[i,j]*dij^3*@f(Cλ,λa[i,j],λr[i,j])*@f(ρ_S)*
    ( (x_0ij^λa[i,j]*(@f(aS_1,λa[i,j])+@f(B,λa[i,j],x_0effij))-
        x_0ij^λr[i,j]*(@f(aS_1,λr[i,j])+@f(B,λr[i,j],x_0effij)))+
        (x_0ij^(λa[i,j]+2)*@f(Q1,λa[i,j])*(@f(aS_1,λa[i,j]+2)+@f(B,λa[i,j]+2,x_0effij))-
         x_0ij^(λr[i,j]+2)*@f(Q1,λr[i,j])*(@f(aS_1,λr[i,j]+2)+@f(B,λr[i,j]+2,x_0effij)))*Dij+
        (x_0ij^(λa[i,j]+4)*@f(Q2,λa[i,j])*(@f(aS_1,λa[i,j]+4)+@f(B,λa[i,j]+4,x_0effij))-
         x_0ij^(λr[i,j]+4)*@f(Q2,λr[i,j])*(@f(aS_1,λr[i,j]+4)+@f(B,λr[i,j]+4,x_0effij)))*Dij^2 )
end
=#

#=
function x_0eff(model::SAFTVRQMieModel, V, T, z, i, j)
    return @f(σeff,i,j)/@f(d,i,j)
end

function a_2(model::SAFTVRQMieModel, V, T, z, i, j,_data = @f(data))
    ϵ = model.params.epsilon.values
    σ = model.params.sigma.values
    λr = model.params.lambda_r.values[i,j]
    λa = model.params.lambda_a.values[i,j]
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)
    x_0ij = @f(x_0,i,j)
    x_0effij = @f(x_0eff,i,j)
    return π*@f(KHS)*(1+@f(χ,i,j))*@f(ρ_S)*ϵ[i,j]^2*@f(d,i,j)^3*@f(C,i,j)^2*(x_0ij^(2*λa)*(@f(aS_1,2*λa)+@f(B,2*λa,x_0effij))-
    x_0ij^(λa+λr)*2*(@f(aS_1,λa+λr)+@f(B,λa+λr,x_0effij))+
    x_0ij^(2*λr)*(@f(aS_1,2*λr)+@f(B,2*λr,x_0effij))+
    x_0ij^(2*λa+2)*2*@f(Q1,λa)*(@f(aS_1,2*λa+2)+@f(B,2*λa+2,x_0effij))*Dij+
    x_0ij^(2*λr+2)*2*@f(Q1,λr)*(@f(aS_1,2*λr+2)+@f(B,2*λr+2,x_0effij))*Dij-
    x_0ij^(λa+λr+2)*2*(@f(Q1,λa)+@f(Q1,λr))*(@f(aS_1,λa+λr+2)+@f(B,λa+λr+2,x_0effij))*Dij+
    x_0ij^(2*λa+4)*@f(Q1,λa)^2*(@f(aS_1,2*λa+4)+@f(B,2*λa+4,x_0effij))*Dij^2+
    x_0ij^(2*λr+4)*@f(Q1,λr)^2*(@f(aS_1,2*λr+4)+@f(B,2*λr+4,x_0effij))*Dij^2-
    x_0ij^(λa+λr+4)*(2*@f(Q1,λa)*@f(Q1,λr))*(@f(aS_1,λa+λr+4)+@f(B,λa+λr+4,x_0effij))*Dij^2+
    x_0ij^(2*λa+4)*2*@f(Q2,λa)*(@f(aS_1,2*λa+4)+@f(B,2*λa+4,x_0effij))*Dij^2+
    x_0ij^(2*λr+4)*2*@f(Q2,λr)*(@f(aS_1,2*λr+4)+@f(B,2*λr+4,x_0effij))*Dij^2-
    x_0ij^(λa+λr+4)*2*(@f(Q2,λa)+@f(Q2,λr))*(@f(aS_1,λa+λr+4)+@f(B,λa+λr+4,x_0effij))*Dij^2+
    x_0ij^(2*λa+6)*2*(@f(Q1,λa)*@f(Q2,λa))*(@f(aS_1,2*λa+6)+@f(B,2*λa+6,x_0effij))*Dij^3+
    x_0ij^(2*λr+6)*2*(@f(Q1,λr)*@f(Q2,λr))*(@f(aS_1,2*λr+6)+@f(B,2*λr+6,x_0effij))*Dij^3-
    x_0ij^(λa+λr+6)*2*(@f(Q1,λr)*@f(Q2,λa)+@f(Q1,λa)*@f(Q2,λr))*(@f(aS_1,λa+λr+6)+@f(B,λa+λr+6,x_0effij))*Dij^3+
    x_0ij^(2*λa+8)*@f(Q2,λa)^2*(@f(aS_1,2*λa+8)+@f(B,2*λa+8,x_0effij))*Dij^4+
    x_0ij^(2*λr+8)*@f(Q2,λr)^2*(@f(aS_1,2*λr+8)+@f(B,2*λr+8,x_0effij))*Dij^4-
    x_0ij^(λa+λr+8)*(2*@f(Q2,λa)*@f(Q2,λr))*(@f(aS_1,λa+λr+8)+@f(B,λa+λr+8,x_0effij))*Dij^4)
end

function χ(model::SAFTVRQMieModel, V, T, z,i,j)
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    σeffij = @f(σeff,i,j)
    ϵeffij = @f(ϵeff,i,j)
    Mwij = (model.params.Mw.values[i] + model.params.Mw.values[j])/2 # check
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)

    ζst_ = ζst(model, V, T, z)
    α = @f(C,i,j)*ϵ[i,j]/ϵeffij*
    (((σ[i,j]/σeffij)^λa[i,j]/(λa[i,j]-3)-(σ[i,j]/σeffij)^λr[i,j]/(λr[i,j]-3))+
        Dij*((σ[i,j]/σeffij)^(2+λa[i,j])*@f(Q1,λa[i,j])/(λa[i,j]-1)-
             (σ[i,j]/σeffij)^(2+λr[i,j])*@f(Q1,λr[i,j])/(λr[i,j]-1))+
        Dij^2*((σ[i,j]/σeffij)^(4+λa[i,j])*@f(Q2,λa[i,j])/(λa[i,j]+1)-
               (σ[i,j]/σeffij)^(4+λr[i,j])*@f(Q2,λr[i,j])/(λr[i,j]+1)))
    return @f(f,α,1)*ζst_+@f(f,α,2)*ζst_^5+@f(f,α,3)*ζst_^8
end



function a_3(model::SAFTVRQMieModel, V, T, z, i, j)
    λr = model.params.lambda_r.values
    λa = model.params.lambda_a.values
    σ = model.params.sigma.values
    ϵ = model.params.epsilon.values
    σeffij = @f(σeff,i,j)
    ϵeffij = @f(ϵeff,i,j)
    Mwij = model.params.Mw.values[i,j]
    Dij = ħ^2/(12*k_B*T*Mwij/N_A*σ[i,j]^2)

    ζst_ = ζst(model, V, T, z)
    α = @f(C,i,j)*ϵ[i,j]/ϵeffij*
    (((σ[i,j]/σeffij)^λa[i,j]/(λa[i,j]-3)-(σ[i,j]/σeffij)^λr[i,j]/(λr[i,j]-3))+
        Dij*((σ[i,j]/σeffij)^(2+λa[i,j])*@f(Q1,λa[i,j])/(λa[i,j]-1)-
             (σ[i,j]/σeffij)^(2+λr[i,j])*@f(Q1,λr[i,j])/(λr[i,j]-1))+
        Dij^2*((σ[i,j]/σeffij)^(4+λa[i,j])*@f(Q2,λa[i,j])/(λa[i,j]+1)-
               (σ[i,j]/σeffij)^(4+λr[i,j])*@f(Q2,λr[i,j])/(λr[i,j]+1)))
    return -ϵeffij^3*@f(f,α,4)*ζst_*exp(@f(f,α,5)*ζst_+@f(f,α,6)*ζst_^2)
end
=#
function a_disp(model::SAFTVRQMieModel, V, T, z, _data = @f(data))
    _σeff, _εff, vrdata = _data
    _d, _ρ_S, ηi, η̅, _η̅st, _, m̄ = vrdata
    nc = length(model)
    Σz    = sum(z)
    m     = model.params.segment.values
    _ε    = model.params.epsilon.values
    _λr   = model.params.lambda_r.values
    _λa   = model.params.lambda_a.values
    _σ    = model.params.sigma.values
    Mw    = model.params.Mw.values
    m̄inv  = 1/m̄
    second_order = (model.fh_order === :fh2)

    a₁ = zero(V + T + first(z));  a₂ = a₁;  a₃ = a₁
    _KHS = SAFTVRMieUtils.KHS(η̅)

    for i in 1:nc
        x_Si = z[i]*m[i]*m̄inv
        for j in 1:i
            x_Sj = z[j]*m[j]*m̄inv

            ε,  λa, λr, σ, dij = _ε[i,j],  _λa[i,j], _λr[i,j], _σ[i,j], _d[i,j]
            σeff, εff          = _σeff[i,j], _εff[i,j]
            C_ij = Cλ_mie(λa, λr)
            dij3 = dij^3
            x_0ij    = σ/dij
            x_0effij = σeff/dij
            Dij = ħ^2/(12*k_B*T*Mw[i,j]/N_A*σ^2)

            # ── a1 : 2π ε dij³ C ρS [ F(λa) − F(λr) ] ─────────────────
            a1_ij = 2π * ε * dij3 * C_ij * _ρ_S *
                    (vrq_F(λa, x_0ij, x_0effij, η̅, Dij, second_order) -
                     vrq_F(λr, x_0ij, x_0effij, η̅, Dij, second_order))

            # ── shared α (used by a2 via χ, and by a3) ────────────────
            Q1λa = vrq_Q1(λa);  Q1λr = vrq_Q1(λr)
            Q2λa = vrq_Q2(λa, second_order)
            Q2λr = vrq_Q2(λr, second_order)
            σcoeff = σ/σeff
            α = C_ij * ε / εff *
                ((σcoeff^λa/(λa-3) - σcoeff^λr/(λr-3)) +
                 Dij   * (σcoeff^(2+λa)*Q1λa/(λa-1) -
                          σcoeff^(2+λr)*Q1λr/(λr-1)) +
                 Dij^2 * (σcoeff^(4+λa)*Q2λa/(λa+1) -
                          σcoeff^(4+λr)*Q2λr/(λr+1)))
            f1, f2, f3, f4, f5, f6 = SAFTVRMieUtils.f123456(α)
            _χ = _χ_from_f123(f1, f2, f3, _η̅st)

            # ── a2 : π KHS (1+χ) ρS ε² dij³ C² (e_λa−e_λr)ᵀ S (e_λa−e_λr) ─
            a2_ij = π * _KHS * (1 + _χ) * _ρ_S * ε^2 * dij3 * C_ij^2 *
                    (vrq_S(λa, λa, x_0ij, x_0effij, η̅, Dij, second_order) -
                     2*vrq_S(λa, λr, x_0ij, x_0effij, η̅, Dij, second_order) +
                     vrq_S(λr, λr, x_0ij, x_0effij, η̅, Dij, second_order))

            # ── a3 ────────────────────────────────────────────────────
            a3_ij = -εff^3 * f4 * _η̅st * exp(_η̅st*(f5 + f6*_η̅st))
            w   = (i == j) ? 1 : 2
            wxy = w * x_Si * x_Sj
            a₁ += a1_ij * wxy
            a₂ += a2_ij * wxy
            a₃ += a3_ij * wxy
        end
    end

    a₁ = a₁ * m̄ / T / Σz
    a₂ = a₂ * m̄ / (T*T) / Σz
    a₃ = a₃ * m̄ / (T*T*T) / Σz
    return a₁ + a₂ + a₃
end

# x_0ij^λ · (aS_1(λ, η̅) + B(λ, x_0effij, η̅))  — η̅ always passed
@inline function vrq_Aexp(λ, x_0ij, x_0effij, η)
    return x_0ij^λ * (SAFTVRMieUtils.aS1(λ, η) + SAFTVRMieUtils.B(λ, x_0effij,η))
end

# F(λ) = Σ_{k=0..2} Dij^k Q_k(λ) x_0ij^(λ+2k) A(λ+2k)          — used by a1
@inline function vrq_F(λ, x_0ij, x_0effij, η̅, Dij, second_order)
    Q1λ = vrq_Q1(λ)
    Q2λ = vrq_Q2(λ, second_order)
    s  = vrq_Aexp(λ,   x_0ij, x_0effij, η̅)
    s += Dij   * Q1λ * vrq_Aexp(λ+2, x_0ij, x_0effij, η̅)
    if second_order
        s += Dij^2 * Q2λ * vrq_Aexp(λ+4, x_0ij, x_0effij, η̅)
    end
    return s
end

# S(λ₁,λ₂) = Σ_{m=0..4} Dij^m C_m(λ₁,λ₂) x_0ij^(λ₁+λ₂+2m) A(λ₁+λ₂+2m)
@inline function vrq_S(λ₁, λ₂, x_0ij, x_0effij, η̅, Dij, second_order)
    Q1₁ = vrq_Q1(λ₁)
    Q1₂ = vrq_Q1(λ₂)
    Q2₁ = vrq_Q2(λ₁, second_order)
    Q2₂ = vrq_Q2(λ₂, second_order)
    λ   = λ₁ + λ₂
    s  = vrq_Aexp(λ,   x_0ij, x_0effij, η̅)
    s += Dij   * (Q1₁ + Q1₂)           * vrq_Aexp(λ+2, x_0ij, x_0effij, η̅)
    s += Dij^2 * (Q2₁ + Q1₁*Q1₂ + Q2₂) * vrq_Aexp(λ+4, x_0ij, x_0effij, η̅)
    if second_order
        s += Dij^3 * (Q1₁*Q2₂ + Q2₁*Q1₂)   * vrq_Aexp(λ+6, x_0ij, x_0effij, η̅)
        s += Dij^4 * (Q2₁*Q2₂)             * vrq_Aexp(λ+8, x_0ij, x_0effij, η̅)
    end
    return s
end

@inline function _χ_from_f123(f1, f2, f3, η)
    η2 = η*η
    η4 = η2*η2
    η5 = η4*η
    η8 = η4*η4
    f1*η + f2*η5 + f3*η8
end

function x0_volume_liquid(model::SAFTVRQMieModel,T,z)
    m = model.params.segment.values
    m̄ = dot(z, m)
    m̄inv = 1/m̄
    comps = @comps
    V = 0.0
    _σ = @f(σeff)
    ηmax = 0.55 #maximum CS
    σ3 = zero(V+T+first(z))
    for i ∈ comps
        x_Si = z[i]*m[i]*m̄inv
        σ3 += x_Si*x_Si*(_σ[i,i]^3)
        for j ∈ 1:i-1
            x_Sj = z[j]*m[j]*m̄inv
            σ3 += 2*x_Si*x_Sj*(_σ[i,j]^3)
        end
    end
    return N_A/ηmax*m̄ * σ3 * π/6
end