struct EmpiricSingleFluidIdealParam <:EoSParam
    a1::Float64
    a2::Float64
    c0::Float64
    u::Vector{Float64}
    v::Vector{Float64}
end

struct EmpiricSingleFluidResidualParam <: EoSParam
    iterators::Vector{UnitRange{Int}}
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Int}
    l::Vector{Int}
    eta::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    epsilon::Vector{Float64}
    b_assoc::Vector{Float64}

    function EmpiricSingleFluidResidualParam(n,t,d,l = Int[],eta = Float64[],beta = Float64[],gamma = Float64[], epsilon = Float64[],b_assoc = Float64[])
        param = new(Vector{UnitRange{Int}}(undef,0),n,t,d,l,eta,beta,gamma,epsilon,b_assoc)
        _calc_iterators!(param)
        return param
    end
end

function _calc_iterators!(param::EmpiricSingleFluidResidualParam)
    n,t,d,l = param.n,param.t,param.d,param.l
    eta,beta,gamma,epsilon = param.eta,param.beta,param.gamma,param.epsilon
    b_assoc = param.b_assoc

    @assert length(n) == length(t) == length(d)
    @assert length(l) < length(d)
    @assert length(eta) == length(beta) == length(gamma) == length(epsilon)
    @assert length(b_assoc) < length(beta)

    #we start from the assoc term, backwards
    length_n = length(n)
    length_beta = length(beta)
    length_b = length(b_assoc)

    length_pol = length_n - length_beta - length(l)
    length_exp = length_n - length_beta
    length_gauss = length_n - length_b
    k_pol = 1:length_pol
    k_exp = (length_pol+1):length_exp
    k_gauss = (length_exp+1):length_gauss
    k_assoc = (length_gauss+1):length_n
    resize!(param.iterators,4)
    param.iterators .= (k_pol,k_exp,k_gauss,k_assoc)
    return param
end

struct EmpiricSingleFluidProperties <: EoSParam
    Mw::Float64 #Molecular Weight, g/mol
    Tc::Float64 #Critical temperature, K
    Pc::Float64 #Critical Pressure,Pa
    rhoc::Float64 #Critical density, mol/m3
    lb_volume::Float64 #lower bound volume, mol/m3
    Ttp::Float64 #triple point temperature, K
    ptp::Float64 #triple point pressure, Pa
    rhov_tp::Float64 #triple point vapor volume, mol/m3
    rhol_tp::Float64 #triple point liquid volume, mol/m3
    acentricfactor::Float64 #acentric factor
    Rgas::Float64 #gas constant used

    function EmpiricSingleFluidProperties(Mw,Tc,Pc,rhoc,lb_volume,
        Ttp = NaN,ptp = NaN, rhov_tp = NaN,rhol_tp = NaN, acentric_factor = NaN, Rgas = R̄)
        return new(Mw,Tc,Pc,rhoc,lb_volume, Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)
    end
end

const ESFProperties = EmpiricSingleFluidProperties
const ESFIdealParam = EmpiricSingleFluidIdealParam
const ESFResidualParam = EmpiricSingleFluidResidualParam


struct EmpiricSingleFluid{𝔾,𝔸} <: EmpiricHelmholtzModel
    type::𝔾
    components::Vector{String}
    properties::ESFProperties
    ancilliaries::𝔸
    ideal::ESFIdealParam
    residual::ESFResidualParam
    references::Vector{String}
end

"""
Single Multiparameter Fluid Equation of state.

```
δ = ρ/ρc
τ = T/Tc
a⁰(δ,τ)   =  log(δ) + a₁ + a₂τ + (c₀ - 1)*log(τ) + ∑vᵢ(1-exp(uᵢτ))
aʳ(δ,τ)   =  aʳ₁+ aʳ₂ + aʳ₃
aʳ₁(δ,τ)  =  ∑nᵢδ^(dᵢ)τ^(tᵢ), i ∈ k_pol
aʳ₂(δ,τ)  =  ∑nᵢexp(-δ^cᵢ)δ^(dᵢ)τ^(tᵢ), i ∈ k_exp
aʳ₃(δ,τ)  =  ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - βᵢ(τ - γᵢ)^2)δ^(dᵢ)τ^(tᵢ), i ∈ k_gauss
aʳ₃(δ,τ)  =  ∑nᵢexp(-ηᵢ(δ - εᵢ)^2 - 1/(βᵢ*(τ -γᵢ)^2 + bᵢ))δ^(dᵢ)τ^(tᵢ), i ∈ k_assoc
```

All parameters are fitted, to allow a equation of state of a single fluid with property calculations as close as possible to the experimental values.
"""
EmpiricSingleFluid

struct IdealEmpiricSingleFluid{𝔾} <: IdealModel
    type::𝔾
    components::Vector{String}
    properties::EmpiricSingleFluidProperties
    ideal::EmpiricSingleFluidIdealParam
    references::Vector{String}
end

function recombine_impl!(model::EmpiricSingleFluid)
    _calc_iterators!(model.residual)
    return model
end

function IdealEmpiricSingleFluid(model::EmpiricSingleFluid)
    return IdealEmpiricSingleFluid(model.type,model.components,model.properties,model.ideal,model.references)
end

idealmodel(model::EmpiricSingleFluid) = IdealEmpiricSingleFluid(model)

R_gas(model::EmpiricSingleFluid) = model.properties.Rgas
R_gas(model::IdealEmpiricSingleFluid) = model.properties.Rgas

function _f0(model::Union{EmpiricSingleFluid,IdealEmpiricSingleFluid},δ,τ)
    a₁ = model.ideal.a1
    a₂ = model.ideal.a2
    c₀ = model.ideal.c0
    α₀ = log(δ) + a₁ + a₂*τ + (c₀ - 1)*log(τ)

    u = model.ideal.u
    length(u) == 0 && return α₀
    v = model.ideal.v
    for i in eachindex(u)
        α₀ += v[i]*log(1 - exp(-u[i]*τ))
    end
    return α₀
end

function _fr1(model::EmpiricSingleFluid,δ,τ,type = model.type)

    αᵣ = zero(δ+τ)
    lnδ = log(δ)
    lnτ = log(τ)

    ℙ = model.residual

    n,t,d,l,η,β,γ,ε,b = ℙ.n,ℙ.t,ℙ.d,ℙ.l,ℙ.eta,ℙ.beta,ℙ.gamma,ℙ.epsilon,ℙ.b_assoc

    k_pol,k_exp,k_gauss,k_assoc = model.residual.iterators

    #strategy for storing.
    #n, t, d, gauss values, always require views
    #l, b does not require views. they are used just once.

    #Polynomial terms
    n_pol = view(n,k_pol)
    t_pol = view(t,k_pol)
    d_pol = view(d,k_pol)
    αᵣ += _fr1_pol(δ,τ,lnδ,lnτ,αᵣ,n_pol,t_pol,d_pol)

    #Exponential terms
    length(k_exp) == 0 && return αᵣ
    n_exp = view(n,k_exp)
    t_exp = view(t,k_exp)
    d_exp = view(d,k_exp)
    αᵣ += _fr1_exp(δ,τ,lnδ,lnτ,αᵣ,n_exp,t_exp,d_exp,l)

    #Gaussian-bell-shaped terms
    length(k_gauss) == 0 && return αᵣ
    n_gauss = view(n,k_gauss)
    t_gauss = view(t,k_gauss)
    d_gauss = view(d,k_gauss)
    αᵣ += _fr1_gauss(δ,τ,lnδ,lnτ,αᵣ,n_gauss,t_gauss,d_gauss,η,β,γ,ε)

    #association terms (new)
    length(k_gao) == 0 && return αᵣ
    lb = length(b)
    lη = length(η)
    k_gao2 = (lη - lb + 1):lη
    n_gao = view(n,k_gao)
    t_gao = view(t,k_gao)
    d_gao = view(d,k_gao)
    η_gao = view(η,k_gao2)
    β_gao = view(β,k_gao2)
    γ_gao = view(γ,k_gao2)
    ε_gao = view(ε,k_gao2)
    b_gao = b
    αᵣ += _fr1_gao(δ,τ,lnδ,lnτ,αᵣ,n_gao,t_gao,d_gao,η_gao,β_gao,γ_gao,ε_gao,b_gao)
    return αᵣ
end

_frx(model::EmpiricSingleFluid{Nothing},δ,τ) = 0.0

function a_ideal(model::IdealEmpiricSingleFluid,V,T,z=SA[1.])
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = Tc/T
    return  _f0(model,δ,τ)
end

a_ideal(model::EmpiricSingleFluid,V,T,z=SA[1.]) = a_ideal(idealmodel(model),V,T,z)

function a_res(model::EmpiricSingleFluid,V,T,z=SA[1.])
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
    τ = Tc/T
    return  _fr1(model,δ,τ) + _frx(model,δ,τ)
end

function eos(model::EmpiricSingleFluid, V, T, z=SA[1.0])
    R = R_gas(model)
    Tc = model.properties.Tc
    rhoc = model.properties.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rhoc
    τ = Tc/T
    return N*R*T*(_f0(model,δ,τ)+_fr1(model,δ,τ) + _frx(model,δ,τ))
end

function eos_res(model::EmpiricSingleFluid,V,T,z=SA[1.0])
    R = R_gas(model)
    Tc = model.consts.Tc
    rhoc = model.consts.rhoc
    N = only(z)
    rho = (N/V)
    δ = rho/rho_c
    τ = Tc/T
    return N*R*T*_fr1(model,δ,τ) + _frx(model,δ,τ)
end

mw(model::EmpiricSingleFluid) = SA[model.properties.Mw]

molecular_weight(model::EmpiricSingleFluid,z = @SVector [1.]) = model.consts.Mw*0.001

T_scale(model::EmpiricSingleFluid,z=SA[1.0]) = model.properties.Tc

p_scale(model::EmpiricSingleFluid,z=SA[1.0]) = model.properties.Pc

lb_volume(model::EmpiricSingleFluid,z=SA[1.0]) = model.properties.lb_volume #finally, an eos model that mentions it max density.

Base.length(::EmpiricSingleFluid) = 1

function Base.show(io::IO,mime::MIME"text/plain",model::EmpiricSingleFluid)
    print(io,"Reference Equation of state for $(model.components[1])")
end

function x0_sat_pure(model::EmpiricSingleFluid,T,z=SA[1.0])
    vv = volume(model.ancilliaries.gas,0.0,T,z)
    vl = volume(model.ancilliaries.liquid,0.0,T,z)
    return (vl,vv)
end

function x0_volume_liquid(model::EmpiricSingleFluid,T,z = SA[1.0])
    volume(model.ancilliaries.liquid,0.0,min(T,model.consts.T_c*one(T)),z)
end

x0_psat(model::EmpiricSingleFluid,T,crit=nothing) = saturation_pressure(model.ancilliaries.saturation,T,SaturationCorrelation())[1]

function x0_saturation_temperature(model::EmpiricSingleFluid,p,z=SA[1.0])
    T = saturation_temperature(model.ancilliaries.saturation,p,SaturationCorrelation())[1]
    vl,vv = x0_sat_pure(model,T)
    return (vl,vv)
end

function crit_pure(model::EmpiricSingleFluid)
    Tc = model.properties.Tc
    Vc = 1/model.properties.rhoc
    Pc = model.properties.Pc

    return (Tc,Pc,Vc)
end

function tryparse_units(val,unit)
    result = try
        #unit_parsed = Unitful.@u_str($unit)
        ThermoState.normalize_units(val*unit)
    catch
        val
    end
    return result
end

function xxxx(path::String)
    _path = only(flattenfilepaths(String[],path))

    json_string = read(_path, String)
    data = JSON3.read(json_string)
    #return data

    #init properties
    info = data[:INFO]
    components = [info[:NAME]]

    eos_data = first(data[:EOS])
    st_data = data[:STATES]
    crit = st_data[:critical]
    rhol_tp_data = st_data[:triple_liquid]
    rhov_tp_data = st_data[:triple_vapor]
    Mw = tryparse_units(get(eos_data,:molar_mass,0.0),get(eos_data,:molar_mass_units,""))
    T_c = tryparse_units(get(crit,:T,NaN),get(crit,:T_units,""))
    P_c = tryparse_units(get(crit,:p,NaN),get(crit,:p_units,""))
    rho_c = tryparse_units(get(crit,:rhomolar,NaN),get(crit,:rhomolar_units,""))

    Ttp = tryparse_units(get(eos_data,:Ttriple,NaN),get(eos_data,:Ttriple_units,""))
    ptp =  tryparse_units(get(rhov_tp_data,:p,NaN),get(rhov_tp_data,:p_units,""))
    if isnan(ptp)
        ptp =  tryparse_units(get(rhol_tp_data,:p,NaN),get(rhov_tl_data,:p_units,""))
    end
    rhol_tp  = tryparse_units(get(rhol_tp_data,:rhomolar,NaN),get(rhol_tp_data,:rhomolar_units,""))
    rhov_tp = tryparse_units(get(rhov_tp_data,:rhomolar,NaN),get(rhov_tp_data,:rhomolar_units,""))
    Rgas = tryparse_units(get(eos_data,:gas_constant,R̄),get(eos_data,:gas_constant_units,""))
    acentric_factor = tryparse_units(get(eos_data,:acentric,NaN),get(eos_data,:acentric_units,""))

    #TODO: in the future, maybe max_density could be in the files?
    lb_volume = 1/(1.25*rhol_tp)
    if isnan(lb_volume)
        lb_volume = 1/(3.25*rho_c)
    end

    properties = EmpiricSingleFluidProperties(Mw,T_c,P_c,rho_c,lb_volume,Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)
    #ideal
    ideal = _parse_ideal(eos_data[:alpha0])
    #residual
    residual = _parse_residual(eos_data[:alphar])
    #ancilliaries
    ancilliaries = _parse_ancilliaries(data[:ANCILLARIES])
    ancilliaries.components[1] = components[1]

    references = [eos_data[:BibTeX_EOS]]

    return EmpiricSingleFluid(nothing,components,properties,ancilliaries,ideal,residual,references)
end

function _parse_ideal(id_data)
    a1 = 0.0
    a2 = 0.0
    c0 = 0.0
    u = Float64[]
    v = Float64[]
    for id_data_i in id_data
        if id_data_i[:type] == "IdealGasHelmholtzLead"
            a1 += id_data_i[:a1]
            a2 += id_data_i[:a2]
        elseif id_data_i[:type] == "IdealGasHelmholtzLogTau"
            c0 += id_data_i[:a] + 1
        elseif id_data_i[:type] == "IdealGasHelmholtzPlanckEinstein"
            append!(v,id_data_i[:n])
            append!(u,id_data_i[:t])
        else
            throw(error("Ideal: $(id_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
        end
    end

    return EmpiricSingleFluidIdealParam(a1,a2,c0,u,v)

end

function _parse_residual(res_data)
   #polynomial y exp terms, we will separate those later
   n = Float64[]
   t = Float64[]
   d = Int[]
   l = Int[]

   #gaussian terms
   n_gauss = Float64[]
   t_gauss = Float64[]
   d_gauss = Int[]
   eta = Float64[]
   beta = Float64[]
   gamma = Float64[]
   epsilon = Float64[]

   #gao association terms
   n_gao = Float64[]
   t_gao = Float64[]
   d_gao = Int[]
   eta_gao = Float64[]
   beta_gao = Float64[]
   gamma_gao = Float64[]
   epsilon_gao = Float64[]
   b_gao = Float64[]

   #non-analytic terms for IAPWS95
   NA_A = Float64[]
   NA_B = Float64[]
   NA_C = Int[]
   NA_D = Float64[]
   NA_a = Float64[]
   NA_b = Float64[]
   NA_beta = Float64[]
   NA_n = Float64[]

   for res_data_i in res_data
       if res_data_i[:type] == "ResidualHelmholtzPower"
           append!(n,res_data_i[:n])
           append!(t,res_data_i[:t])
           append!(d,res_data_i[:d])
           append!(l,res_data_i[:l])
       elseif res_data_i[:type] == "ResidualHelmholtzGaussian"
           append!(n_gauss,res_data_i[:n])
           append!(t_gauss,res_data_i[:t])
           append!(d_gauss,res_data_i[:d])
           append!(eta,res_data_i[:eta])
           append!(beta,res_data_i[:beta])
           append!(gamma,res_data_i[:gamma])
           append!(epsilon,res_data_i[:epsilon])
       elseif res_data_i[:type] == "ResidualHelmholtzGaoB"
           append!(n_gao,res_data_i[:n])
           append!(t_gao,res_data_i[:t])
           append!(d_gao,res_data_i[:d])
           append!(eta_gao,res_data_i[:eta])
           append!(beta_gao,res_data_i[:beta])
           append!(gamma_gao,res_data_i[:gamma])
           append!(epsilon_gao,res_data_i[:epsilon])
           append!(b_gao,res_data_i[:b])
       elseif res_data_i[:type] == "ResidualHelmholtzNonAnalytic"
           append!(NA_A,res_data_i[:A])
           append!(NA_B,res_data_i[:B])
           append!(NA_C,res_data_i[:C])
           append!(NA_D,res_data_i[:D])
           append!(NA_a,res_data_i[:a])
           append!(NA_b,res_data_i[:b])
           append!(NA_beta,res_data_i[:beta])
           append!(NA_n,res_data_i[:n])
       else
           throw(error("Residual: $(res_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
       end
   end

   pol_vals = findall(iszero,l)
   exp_vals = findall(!iszero,l)
   _n = vcat(n[pol_vals],n[exp_vals],n_gauss,n_gao)
   _t = vcat(t[pol_vals],t[exp_vals],t_gauss,t_gao)
   _d = vcat(d[pol_vals],d[exp_vals],d_gauss,d_gao)
   _l = l[exp_vals]
   _η = vcat(eta,eta_gao)
   _β = vcat(beta,beta_gao)
   _γ = vcat(gamma,gamma_gao)
   _ε = vcat(gamma,gamma_gao)
   _b = b_gao
   return EmpiricSingleFluidResidualParam(_n,_t,_d,_l,_η,_β,_γ,_ε,_b)
end

function _parse_ancilliaries(anc_data)
    #saturation pressure
    p_data = anc_data[:pS]
    rhol_data = anc_data[:rhoL]
    rhov_data = anc_data[:rhoV]

    ps_anc = if p_data[:type] == "pV"
        T_c = p_data[:T_r]
        P_c = p_data[:reducing_value] * 1.0
        n = p_data[:n]
        t = p_data[:t]
        PolExpSat(T_c,P_c,n,t)
    else
        throw(error("Ancilliary saturation pressure: $(id_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
    end

    rhov_anc = if rhov_data[:type] == "rhoV"
        T_c = rhov_data[:T_r]
        rho_c = rhov_data[:reducing_value] * 1.0
        n = rhov_data[:n]
        t = rhov_data[:t]
        PolExpVapour(T_c,rho_c,n,t)
    else
        throw(error("Ancilliary vapour density: $(id_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
    end

    rhol_anc = if rhol_data[:type] == "rhoLnoexp"
        T_c = rhol_data[:T_r]
        rho_c = rhol_data[:reducing_value] * 1.0
        n = rhol_data[:n]
        t = rhol_data[:t]
        PolExpVapour(T_c,rho_c,n,t)
    else
        throw(error("Ancilliary liquid density: $(id_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
    end

    return CompositeModel(["ancilliaries"],gas = rhov_anc,liquid = rhol_anc,saturation = ps_anc)
end
export EmpiricSingleFluid

