JSON_ALTERNATIVE_NAMES = Dict{String,String}(
    "carbon dioxide" => "CarbonDioxide",
    "hydrogen chloride" => "HydrogenChloride",
    "hydrogen sulfide" => "HydrogenSulfide",
    "isopentane" => "ISOPENTANE",
    "nonane" => "n-Nonane",
    "octane" => "n-Octane",
    "heptane" => "n-Heptane",
    "hexane" => "n-Hexane",
    "pentane" => "n-Pentane",
    "butane" => "n-Butane",
    "propane" => "n-Propane",
    "decane" => "n-Decane",
    "undecane" => "n-Undecane",
    "dodecane" => "n-Dodecane",
    "carbonmonoxide" => "CARBONMONOXIDE",
    "hydrogensulfide" => "HydrogenSulfide",
)

function is_coolprop_loaded()
    lib_handler = Base.Libc.Libdl.dlopen("libcoolprop";throw_error = false)
    !isnothing(lib_handler)
end

function coolprop_csv(component::String)
    lib_handler = Base.Libc.Libdl.dlopen("libcoolprop";throw_error = false)
    if !isnothing(lib_handler)
       #libcoolprop is present.
        buffer_length = 2<<12
        message_buffer = Vector{UInt8}(undef,buffer_length)
        method_handler = Base.Libc.Libdl.dlsym(lib_handler,:get_fluid_param_string)
        err_handler = Base.Libc.Libdl.dlsym(lib_handler,:get_global_param_string)
        val = 0
        for i in 1:5
            val = ccall(method_handler, Clong, (Cstring, Cstring, Ptr{UInt8}, Int), component, "JSON", message_buffer::Array{UInt8, 1}, buffer_length)
            if val == 0
                ccall(err_handler, Clong, (Cstring, Ptr{UInt8}, Int), "errstring", message_buffer::Array{UInt8, 1}, buffer_length)
                err = unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
                if err == "Buffer size is too small"
                    resize!(message_buffer,buffer_length<<1)
                    buffer_length = length(message_buffer)
                else
                    return false,unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
                end
            else
                return true,unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
            end
        end
        return false,unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
    else
        return false,""
    end

end
function tryparse_units(val,unit)
    result = try
        unit_parsed = Unitful.uparse(unit)
        ThermoState.normalize_units(val*unit_parsed)
    catch
        val
    end
    return result
end

get_only_comp(x::Vector{String}) = only(x)
get_only_comp(x::String) = x

function get_json_data(components;
    userlocations = String[],
    coolprop_userlocations = true,
    verbose = false,
    )

    component = get_only_comp(components)
    if first(component) != '{' #not json
        _paths = flattenfilepaths(["Empiric"],userlocations)
        
        norm_comp1 = normalisestring(component)
        f0 = x -> normalisestring(last(splitdir(first(splitext(x))))) == norm_comp1

        found_paths = filter(f0,_paths)
        if iszero(length(found_paths))
            #try to extract from coolprop.
            !coolprop_userlocations && throw(error("cannot found component file $(component)."))
            !is_coolprop_loaded() && throw(error("cannot found component file $(component). Try loading the CoolProp library by loading it."))
            
            alternative_comp = get(JSON_ALTERNATIVE_NAMES,norm_comp1,norm_comp1)
            success,json_string = coolprop_csv(alternative_comp)
            if success
                data = JSON3.read(json_string)[1]
                return data
            else
                if length(json_string) == 0
                    throw(error("cannot found component file $(component)."))
                else
                    throw(error("Coolprop: $(json_string)."))
                end
            end
        end
        _path = last(found_paths)
        json_string = read(_path, String)
        data = JSON3.read(json_string)
    else
        data = JSON3.read(component)
    end
    return data
end

function SingleFluid(components;
        userlocations = String[],
        ancillaries = nothing,
        ancillaries_userlocations = String[],
        estimate_pure = false,
        coolprop_userlocations = true,
        Rgas = nothing,
        verbose = false)


    components = [get_only_comp(components)]
    data = try
        get_json_data(components;userlocations,coolprop_userlocations,verbose)
        catch e
            !estimate_pure && rethrow(e)
            nothing
        end
    if data === nothing && estimate pure
        return XiangDeiters(components;userlocations,verbose = verbose)
    end
    eos_data = first(data[:EOS])
    #properties
    properties = _parse_properties(data,Rgas,verbose)
    #ideal
    ideal = _parse_ideal(eos_data[:alpha0],verbose)
    #residual
    residual = _parse_residual(eos_data[:alphar],verbose)
    #ancillaries
    if ancillaries === nothing
        init_ancillaries = _parse_ancillaries(data[:ANCILLARIES],verbose)
        init_ancillaries.components[1] = components[1]
    else
        init_ancillaries = init_model(ancillaries,components,ancillaries_userlocations,verbose)
    end

    references = [eos_data[:BibTeX_EOS]]

    return SingleFluid(components,properties,init_ancillaries,ideal,residual,references)
end

function SingleFluidIdeal(components;
    userlocations = String[],
    Rgas = nothing,
    verbose = false,
    coolprop_userlocations = true)

    data = get_json_data(components;userlocations,coolprop_userlocations,verbose)
    components = [get_only_comp(components)]
    eos_data = first(data[:EOS])
    #properties
    properties = _parse_properties(data,Rgas,verbose)
    #ideal
    ideal = _parse_ideal(eos_data[:alpha0],verbose)

    references = [eos_data[:BibTeX_EOS]]

    return SingleFluidIdeal(components,properties,ideal,references)
end


function _parse_properties(data,Rgas0 = nothing, verbose = false)
    info = data[:INFO]
    eos_data = first(data[:EOS])
    st_data = data[:STATES]
    crit = st_data[:critical]
    reducing = eos_data[:STATES][:reducing]

    Tr = tryparse_units(get(reducing,:T,NaN),get(reducing,:T_units,""))
    rhor = tryparse_units(get(reducing,:rhomolar,NaN),get(reducing,:rhomolar_units,""))
    Mw = 1000*tryparse_units(get(eos_data,:molar_mass,0.0),get(eos_data,:molar_mass_units,""))
    T_c = tryparse_units(get(crit,:T,NaN),get(crit,:T_units,""))
    P_c = tryparse_units(get(crit,:p,NaN),get(crit,:p_units,""))
    rho_c = tryparse_units(get(crit,:rhomolar,NaN),get(crit,:rhomolar_units,""))

    rhov_tp_data = get(st_data,:triple_vapor,nothing)
    Ttp = tryparse_units(get(eos_data,:Ttriple,NaN),get(eos_data,:Ttriple_units,""))
    if rhov_tp_data !== nothing
        ptp =  tryparse_units(get(rhov_tp_data,:p,NaN),get(rhov_tp_data,:p_units,""))
        rhov_tp = tryparse_units(get(rhov_tp_data,:rhomolar,NaN),get(rhov_tp_data,:rhomolar_units,""))
    else
        ptp,rhov_tp = NaN,NaN
    end

    rhol_tp_data = get(st_data,:triple_liquid,nothing)
    if rhol_tp_data !== nothing
        if isnan(ptp)
            ptp =  tryparse_units(get(rhol_tp_data,:p,NaN),get(rhov_tl_data,:p_units,""))
        end
        rhol_tp  = tryparse_units(get(rhol_tp_data,:rhomolar,NaN),get(rhol_tp_data,:rhomolar_units,""))
    else
        rhol_tp = NaN
    end
    if  Rgas0 === nothing
        Rgas = tryparse_units(get(eos_data,:gas_constant,R̄),get(eos_data,:gas_constant_units,""))
    else
        Rgas = Rgas0
    end
    acentric_factor = tryparse_units(get(eos_data,:acentric,NaN),get(eos_data,:acentric_units,""))

    #TODO: in the future, maybe max_density could be in the files?

    lb_volume = 1/tryparse_units(get(crit,:rhomolar_max,NaN),get(crit,:rhomolar_max_units,""))
    isnan(lb_volume) && (lb_volume = 1/tryparse_units(get(eos_data,:rhomolar_max,NaN),get(eos_data,:rhomolar_max_units,"")))
    isnan(lb_volume) && (lb_volume = 1/(1.25*rhol_tp))
    isnan(lb_volume) && (lb_volume = 1/(3.25*rho_c))
    return SingleFluidProperties(Mw,Tr,rhor,lb_volume,T_c,P_c,rho_c,Ttp,ptp,rhov_tp,rhol_tp,acentric_factor,Rgas)
end

function _parse_ideal(id_data,verbose = false)
    a1 = 0.0
    a2 = 0.0
    c0 = 0.0
    R0 = 0.0
    n = Float64[]
    t = Float64[]
    c = Float64[]
    d = Float64[]
    np = Float64[]
    tp = Float64[]
    n_gerg = Float64[]
    v_gerg = Float64[]
    for id_data_i in id_data
        if id_data_i[:type] == "IdealGasHelmholtzLead" || id_data_i[:type] == "IdealGasHelmholtzEnthalpyEntropyOffset"
            a1 += id_data_i[:a1]
            a2 += id_data_i[:a2]
        elseif id_data_i[:type] == "IdealGasHelmholtzLogTau"
            c0 += id_data_i[:a]
        elseif id_data_i[:type] == "IdealGasHelmholtzPlanckEinstein"
            append!(n,id_data_i[:n])
            append!(t,-id_data_i[:t])
            l = length(id_data_i[:n])
            append!(c,fill(1.,l))
            append!(d,fill(-1.,l))
        elseif id_data_i[:type] == "IdealGasHelmholtzPlanckEinsteinFunctionT"
            _Tc = id_data_i[:Tcrit]
            append!(n,id_data_i[:n])
            append!(t, -id_data_i[:v] ./ _Tc)
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
            a2 += -cpi*τ0
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
            t_pj = id_data_i[:t]
            n_pj = id_data_i[:n]
            for i in 1:length(t)
                #workaround 1: it seems that sometinmes, people store lead as power
                #it is more efficient if we transform from power to lead term, if possible
                if t_pj[i] == 0
                    a1 += n_pj[i]
                elseif t_pj[i] == 1
                    a2 += n_pj[i]
                else
                    push!(np,n_pj[i])
                    push!(tp,t_pj[i])
                end
            end
        elseif id_data_i[:type] == "IdealHelmholtzPlanckEinsteinGeneralized" || id_data_i[:type] == "IdealGasHelmholtzPlanckEinsteinGeneralized"
            append!(n,id_data_i[:n])
            append!(t,id_data_i[:t])
            append!(c,id_data_i[:c])
            append!(d,id_data_i[:d])
        elseif id_data_i[:type] == "IdealGasHelmholtzCP0AlyLee"
            alylee_data = id_data_i[:c]
            _Tc = id_data_i[:Tc]
            _T0 = id_data_i[:T0]

            @assert length(alylee_data) == 5 "aly-lee is defined with only 5 terms. add an additional ally lee term if you require more coefficients."
            A,B,C,D,E = alylee_data

            if !iszero(A)
                τ0 = _T0/_Tc
                a1 += A*(1 - log(τ0))
                a2 += -A*_T0/_Tc
                c0 += A
            end

            push!(n,B)
            push!(t,-2*C/Tc)
            push!(c,1)
            push!(d,-1)

            push!(n,-D)
            push!(t,-2*E/Tc)
            push!(c,1)
            push!(d,1)
        elseif id_data_i[:type] == "IdealGasClapeyronJLGerg2008"
            append!(n_gerg,id_data_i[:n])
            append!(v_gerg,id_data_i[:v])
        elseif id_data_i[:type] == "IdealGasClapeyronJLR0"
            R0 = id_data_i[:R0]
        else
            throw(error("Ideal: $(id_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
        end
    end

    return SingleFluidIdealParam(a1,a2,c0,n,t,c,d,np,tp,n_gerg,v_gerg,R0)

end

function _parse_residual(res_data, verbose = false)
    #polynomial y exp terms, we will separate those later
    n = Float64[]
    t = Float64[]
    d = Int[]
    l = Int[]
    g = Float64[]

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
            append!(g,ones(length(res_data_i[:l])))
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
            append!(n,res_data_i[:n])
            append!(t,res_data_i[:t])
            append!(d,res_data_i[:d])
            append!(l,res_data_i[:l])
            append!(g,res_data_i[:g])
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
    _g = g[exp_vals]
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

   return SingleFluidResidualParam(_n,_t,_d,_l,_g,_η,_β,_γ,_ε;gao_b,na,assoc)
end

function _parse_ancilliary_func(anc,input_key,output_key)
    anc_typemap = Dict(
    "pV" => :exp,
    "pL" => :exp,
    "rhoV" => :exp,
    "rhoL" => :exp,
    "rhoLnoexp" => :noexp,
    "rhoVnoexp" => :noexp,
    "rational" => :rational,
    )

    anc_using_r_map = Dict(
    "pV" => true,
    "pL" => true,
    "rhoV" => false,
    "rhoL" => false,
    "rhoLnoexp" => false,
    "rhoVnoexp" => false,
    "rational" => false,
    )

    input_r = anc[input_key] * 1.0
    output_r = anc[output_key] * 1.0
    n = Float64.(anc[:n])
    t = Float64.(anc[:t])
    type_str = anc[:type]
    using_input_r = get(anc,:using_tau_r,anc_using_r_map[type_str])
    type = get(anc_typemap,anc[:type],Symbol(type_str))
    return GenericAncEvaluator(n,t,input_r,output_r,type,using_input_r)
end

function _parse_ancillaries(anc_data,verbose = false)
    #saturation pressure
    p_data = anc_data[:pS]
    rhol_data = anc_data[:rhoL]
    rhov_data = anc_data[:rhoV]

    ps_anc = PolExpSat(_parse_ancilliary_func(p_data,:T_r,:reducing_value))
    rhov_anc = PolExpVapour(_parse_ancilliary_func(rhov_data,:T_r,:reducing_value))
    rhol_anc = PolExpVapour(_parse_ancilliary_func(rhol_data,:T_r,:reducing_value))
    return CompositeModel(["ancillaries"],gas = rhov_anc,liquid = rhol_anc,saturation = ps_anc)
end
export SingleFluid

function allxxx()
    _path = flattenfilepaths("Empiric/test",String[])
    res = String[]

    for p in _path
        json_string = read(p, String)
        data = JSON3.read(json_string)
        a0data = data[:ANCILLARIES]
        if haskey(a0data,:pS)
            !a0data[:pS][:using_tau_r] && println(data[:INFO][:NAME])
            push!(res,a0data[:pS][:type])
        else
            #pV: "p'' = pc*exp(Tc/T*sum(n_i*theta^t_i))"
            #
        end
        #for a0 in a0data
        #    a0type = a0[:type]
        #    #a0type == "ResidualHelmholtzNonAnalytic" && println(data[:INFO][:NAME])
        #    push!(res,a0[:type])
        #end
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
 `IdealGasHelmholtzCP0AlyLee` done, needed for GERG2008, but instead our own type is used.
 `IdealGasHelmholtzCP0Constant` done

all residual types

 `ResidualHelmholtzPower` done
 `ResidualHelmholtzGaussian` done
 `ResidualHelmholtzGaoB` done
 `ResidualHelmholtzNonAnalytic` done, only water have it, maybe optimize the heck out of it?
 `ResidualHelmholtzExponential` not done: (Fluorine,Propyne,R114,R13,R14,R21,RC318)
 `ResidualHelmholtzAssociating` not done, only methanol have it
 `ResidualHelmholtzLemmon2005` mpt done, only R125 have it
=#

#=
AlyLee parser:

[A,B,C,D,E]

push!(n_gpe,B)
push!(t_gpe,-2*C/Tc)
push!(c_gpe,1)
push!(d_gpe,-1)

push!(n_gpe,-D)
push!(t_gpe,-2*E/Tc)
push!(c_gpe,1)
push!(d_gpe,1)
=#