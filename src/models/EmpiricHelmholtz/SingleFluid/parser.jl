function tryparse_units(val,unit)
    result = try
        unit_parsed = Unitful.uparse(unit)
        ThermoState.normalize_units(val*unit_parsed)
    catch
        val
    end
    return result
end

function fff(path::String)
    _path = only(flattenfilepaths(String[],path))

    json_string = read(_path, String)
    data = JSON3.read(json_string)
end

function SingleFluid(component::String;userlocations = String[])
    _paths = flattenfilepaths(["Empiric","Empiric/test"],userlocations)
    normalized_comp = normalisestring(component)
    f0 = x -> normalisestring(last(splitdir(first(splitext(x))))) == normalized_comp
    _path = last(filter(f0,_paths))
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
    Mw = 1000*tryparse_units(get(eos_data,:molar_mass,0.0),get(eos_data,:molar_mass_units,""))
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
    
    lb_volume = tryparse_units(get(crit,:rhomolar_max,NaN),get(crit,:rhomolar_max_units,""))
    isnan(lb_volume) && (lb_volume = 1/(1.25*rhol_tp))
    isnan(lb_volume) && (lb_volume = 1/(3.25*rho_c))

    properties = EmpiricSingleFluidProperties(Mw,T_c,P_c,rho_c,lb_volume,Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)
    #ideal
    ideal = _parse_ideal(eos_data[:alpha0])
    #residual
    residual = _parse_residual(eos_data[:alphar])
    #ancilliaries
    ancilliaries = _parse_ancilliaries(data[:ANCILLARIES])
    ancilliaries.components[1] = components[1]

    references = [eos_data[:BibTeX_EOS]]

    return EmpiricSingleFluid(components,properties,ancilliaries,ideal,residual,references)
end

function _parse_ideal(id_data)
    a1 = 0.0
    a2 = 0.0
    c0 = 0.0
    n = Float64[]
    t = Float64[]
    c = Float64[]
    d = Float64[]
    np = Float64[]
    tp = Float64[]
    for id_data_i in id_data
        if id_data_i[:type] == "IdealGasHelmholtzLead" || id_data_i[:type] == "IdealGasHelmholtzEnthalpyEntropyOffset"
            a1 += id_data_i[:a1]
            a2 += id_data_i[:a2]
        elseif id_data_i[:type] == "IdealGasHelmholtzLogTau"
            c0 += id_data_i[:a]
        elseif id_data_i[:type] == "IdealGasHelmholtzPlanckEinstein"
            append!(n,id_data_i[:n])
            append!(t,id_data_i[:t])
            l = length(id_data_i[:n])
            append!(c,fill(1.,l))
            append!(d,fill(-1.,l))
        elseif id_data_i[:type] == "IdealGasHelmholtzCP0Constant"
            #c - cT0/Tc*τ - c*(log(τ/τ0))
            #c - cT0/Tc*τ - c*(log(τ) - log(τ0))
            #c + c*log(τ0) - cT0/Tc*τ - c*(log(τ))
            cpi = id_data_i[:cp_over_R]
            _T0 = id_data_i[:T0]
            _Tc = id_data_i[:Tc]
            τ0 = _T0/_Tc
            a1 += cpi*(1 - log(τ0))
            a2 += -cpi*_T0/_Tc
            c0 += cpi
        elseif id_data_i[:type] == "IdealGasHelmholtzCP0PolyT"
            _T0 = id_data_i[:T0]
            _Tc = id_data_i[:Tc]
            cp_c = id_data_i[:c]
            cp_t = id_data_i[:t]
            τ0 = _T0/_Tc

            for i in eachindex(cp_t)
                ti = t[i]
                cpi = cp_c[i]
                if ti == 0
                    a1 += cpi*(1 - log(τ0))
                    a2 += -cpi*_T0/_Tc
                    c0 += cpi
                else
                    T0t = _T0^ti
                    a1 += (cpi*T0t)/t
                    a2 += (-T0t*_T0)*cpi/(_Tc*(ti+1))
                    push!(tp,-ti)
                    push!(np,(c*(_Tc^ti))*(1/(ti+1) - 1/ti))
                end
            end
        elseif id_data_i[:type] == "IdealGasHelmholtzPower"
            t = id_data_i[:t]
            n = id_data_i[:n]
            for i in 1:length(t)
                #workaround 1: it seems that sometinmes, people store lead as power
                #it is more efficient if we transform from power to lead term, if possible
                if t[i] == 0
                    a1 += n[i]
                elseif t[i] == 1
                    a2 += n[i]
                else
                    push!(np,n[i])
                    push!(tp,t[i])
                end
            end
        elseif id_data_i[:type] == "IdealHelmholtzPlanckEinsteinGeneralized"
            append!(n,id_data_i[:n])
            append!(t,id_data_i[:t])
            append!(c,id_data_i[:c])
            append!(d,id_data_i[:d])
        else
            throw(error("Ideal: $(id_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
        end
    end

    return EmpiricSingleFluidIdealParam(a1,a2,c0,n,t,c,d,np,tp)

end

function _parse_residual(res_data)
    #polynomial y exp terms, we will separate those later
    n = Float64[]
    t = Float64[]
    d = Int[]
    l = Int[]

    #modified exponential terms
    exp_n = Float64[]
    exp_t = Float64[]
    exp_d = Float64[]
    exp_l = Float64[]
    exp_gamma = Float64[]

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

    #assoc terms
    assoc_epsilonbar = 0.0
    assoc_kappabar = 0.0
    assoc_a = 0.0 
    assoc_m = 0.0
    assoc_vbarn = 0.0
    assoc = false

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
        elseif res_data_i[:type] == "ResidualHelmholtzExponential"
            append!(exp_n,res_data_i[:n])
            append!(exp_t,res_data_i[:t])
            append!(exp_d,res_data_i[:d])
            append!(exp_l,res_data_i[:l])
            append!(exp_gamma,res_data_i[:gamma])
        elseif res_data_i[:type] == "ResidualHelmholtzAssociating"
            if assoc == true
                throw(error("Residual: $(res_data_i[:type]) we only support one Associating term."))
            end
            assoc = true
            assoc_epsilonbar += res_data_i[:epsilonbar]
            assoc_kappabar += res_data_i[:kappabar]
            assoc_a += res_data_i[:a]
            assoc_m += res_data_i[:m]
            assoc_vbarn += res_data_i[:vbarn]
        else
            throw(error("Residual: $(res_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
        end
    end

    pol_vals = findall(iszero,l)
    exp_vals = findall(!iszero,l)
    _n = vcat(n[pol_vals],n[exp_vals],n_gauss)
    _t = vcat(t[pol_vals],t[exp_vals],t_gauss)
    _d = vcat(d[pol_vals],d[exp_vals],d_gauss)
    _l = l[exp_vals]
    _η = eta
    _β = beta
    _γ = gamma
    _ε = epsilon
    
    #gao_b term
    gao_b = GaoBTerm(n_gao,t_gao,d_gao,eta_gao,beta_gao,gamma_gao,epsilon_gao,b_gao)

    #non analytical term
    na = NonAnalyticTerm(NA_A,NA_B,NA_C,NA_D,NA_a,NA_b,NA_beta,NA_n)
    
    #assoc terms
    assoc = Associating2BTerm(assoc_epsilonbar,assoc_kappabar,assoc_a,assoc_m,assoc_vbarn)
    
    #exponential term
    exp = ExponentialTerm(exp_n,exp_t,exp_d,exp_l,exp_gamma)

   return EmpiricSingleFluidResidualParam(_n,_t,_d,_l,_η,_β,_γ,_ε;gao_b,na,assoc,exp)
end

function _parse_ancilliaries(anc_data)
    #saturation pressure
    p_data = anc_data[:pS]
    rhol_data = anc_data[:rhoL]
    rhov_data = anc_data[:rhoV]

    ps_anc = if p_data[:type] in ("pV","pL")
        T_c = p_data[:T_r]
        P_c = p_data[:reducing_value] * 1.0
        n = Float64.(p_data[:n])
        t = Float64.(p_data[:t])
        PolExpSat(T_c,P_c,n,t)
    else
        throw(error("Ancilliary saturation pressure: $(p_data[:type]) not supported for the moment. open an issue in the repository for help."))
    end

    rhov_anc = if rhov_data[:type] == "rhoV"
        T_c = rhov_data[:T_r]
        rho_c = rhov_data[:reducing_value] * 1.0
        n = Float64.(p_data[:n])
        t = Float64.(p_data[:t])
        PolExpVapour(T_c,rho_c,n,t)
    else
        throw(error("Ancilliary vapour density: $(rhov_data[:type]) not supported for the moment. open an issue in the repository for help."))
    end

    rhol_anc = if rhol_data[:type] == "rhoLnoexp"
        T_c = rhol_data[:T_r]
        rho_c = rhol_data[:reducing_value] * 1.0
        n = Float64.(p_data[:n])
        t = Float64.(p_data[:t])
        PolExpVapour(T_c,rho_c,n,t)
    else
        throw(error("Ancilliary liquid density: $(rhol_data[:type]) not supported for the moment. open an issue in the repository for help."))
    end

    return CompositeModel(["ancilliaries"],gas = rhov_anc,liquid = rhol_anc,saturation = ps_anc)
end
export EmpiricSingleFluid

function allxxx()
    _path = flattenfilepaths("Empiric/test",String[])
    res = String[]
    ct = 0
    for p in _path
        json_string = read(p, String)
        data = JSON3.read(json_string)
        a0data = data[:EOS][1][:alphar]
        for a0 in a0data
            a0type = a0[:type]
            a0type == "ResidualHelmholtzNonAnalytic" && println(data[:INFO][:NAME])
            push!(res,a0[:type])
        end
    end
    return unique!(res)
end

#=

all ideal types

 `IdealGasHelmholtzLead` done
 `IdealGasHelmholtzLogTau` done
 `IdealGasHelmholtzPlanckEinstein` done
 `IdealGasHelmholtzEnthalpyEntropyOffset` done, same as Lead
 `IdealGasHelmholtzPower` done
 `IdealGasHelmholtzPlanckEinsteinGeneralized` done
 `IdealGasHelmholtzCP0PolyT` done
 `IdealGasHelmholtzCP0AlyLee` #not done, only n-Heptane and D6 have it, only 5 terms
 `IdealGasHelmholtzCP0Constant` done

all residual types

 `ResidualHelmholtzPower` done
 `ResidualHelmholtzGaussian` done
 `ResidualHelmholtzGaoB` done
 `ResidualHelmholtzNonAnalytic` done, oly water have it, maybe optimize the heck out of it?
 `ResidualHelmholtzExponential` not done: (Fluorine,Propyne,R114,R13,R14,R21,RC318)
 `ResidualHelmholtzAssociating` not done, only methanol have it
 `ResidualHelmholtzLemmon2005` mpt done, only R125 have it
=#