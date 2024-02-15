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

function coolprop_crit_data end

@static if !isdefined(Base,:get_extension)
    function coolprop_handler()
        #for some reason, this does not work on linux/mac
        lib_handler1 = Base.Libc.Libdl.dlopen(:libcoolprop;throw_error = false)
        #return lib_handler1
        lib_handler1 !== nothing && return lib_handler1
        if !Sys.iswindows()
            #search on all dynamic libs, filter libCoolProp. TODO: find something faster.
            dllist = Base.Libc.Libdl.dllist()
            x =findall(z->occursin("libCoolProp",z),dllist)
            length(x) == 0 && return nothing
            t = dllist[x[1]]
            lib_handler2 = Base.Libc.Libdl.dlopen(t;throw_error = false)
            return lib_handler2
        else
            return lib_handler1
        end
    end
else
    #defined in ClapeyronCoolPropExt
    function coolprop_handler end
end

function is_coolprop_loaded()
    handler = coolprop_handler()
    res = handler !== nothing
    Base.Libc.Libdl.dlclose(handler)
    return res
end

function coolprop_csv(component::String,comp = "")
    lib_handler = coolprop_handler()
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
                    Base.Libc.Libdl.dlclose(lib_handler)
                    return false,unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
                end
            else
                Base.Libc.Libdl.dlclose(lib_handler)
                return true,unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
            end
        end
        Base.Libc.Libdl.dlclose(lib_handler)
        return false,unsafe_string(convert(Ptr{UInt8}, pointer(message_buffer::Array{UInt8, 1})))
    else
        Base.Libc.Libdl.dlclose(lib_handler)
        throw(error("cannot found component file $(comp). Try loading the CoolProp library by loading it."))
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
            verbose && @info "JSON for $(info_color(component)) not found in supplied paths"
            verbose && coolprop_userlocations && @info "trying to look JSON for $(info_color(component)) in CoolProp"
            #try to extract from coolprop.
            !coolprop_userlocations && throw(error("cannot found component file $(component)."))
            alternative_comp = get(JSON_ALTERNATIVE_NAMES,norm_comp1,component)
            success,json_string = coolprop_csv(alternative_comp,component)
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
        verbose && @info "JSON found: $_path"

        json_string = read(_path, String)
        data = JSON3.read(json_string)
    else
        verbose && @info "parsing supplied JSON data."
        data = JSON3.read(component)
    end
    return data
end

"""
    SingleFluid(components;
            userlocations = String[],
            ancillaries = nothing,
            ancillaries_userlocations = String[],
            estimate_pure = false,
            coolprop_userlocations = true,
            Rgas = nothing,
            verbose = false)

## Input parameters
- JSON data (CoolProp and teqp format)

## Input models
- `ancillaries`: a model that provides initial guesses for saturation calculations. if `nothing`, then they will be parsed from the input JSON.

## Description

Instantiates a single-component Empiric EoS model. `Rgas` can be used to set the value of the gas constant that is used during property calculations.

If `coolprop_userlocations` is true, then Clapeyron will try to look if the fluid is present in the CoolProp library.

The properties, ideal and residual terms can be accessed via the `properties`, `ideal` and `residual` fields respectively:

```julia-repl
julia> model = SingleFluid("water")
MultiParameter Equation of state for water:
 Polynomial power terms: 7
 Exponential terms: 44
 Gaussian bell-shaped terms: 3
 Non Analytic terms: 2

julia> model.ideal
Ideal MultiParameter coefficients:
 Lead terms: -8.3204464837497 + 6.6832105275932*τ + 3.00632*log(τ)
 Plank-Einstein terms: 5

julia> model.residual
Residual MultiParameter coefficients:
 Polynomial power terms: 7
 Exponential terms: 44
 Gaussian bell-shaped terms: 3
 Non Analytic terms: 2
```
"""
SingleFluid

function SingleFluid(components;
        userlocations = String[],
        ancillaries = nothing,
        ancillaries_userlocations = String[],
        estimate_pure = false,
        coolprop_userlocations = true,
        Rgas = nothing,
        verbose = false,
        idealmodel = nothing,
        ideal_userlocations = String[])


    _components = format_components(components)
    single_component_check(SingleFluid,_components)

    data = try
        get_json_data(_components;userlocations,coolprop_userlocations,verbose)
        catch e
            !estimate_pure && rethrow(e)
            nothing
        end
    if data === nothing && estimate_pure
        return XiangDeiters(components;userlocations,verbose = verbose,idealmodel = idealmodel,ideal_userlocations = ideal_userlocations)
    end
    eos_data = first(data[:EOS])
    #properties
    properties = _parse_properties(data,Rgas,verbose)
    #ideal
    if idealmodel === nothing
        ideal_data = eos_data[:alpha0]
    else
        init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
        ideal_data = Clapeyron.idealmodel_to_json_data(init_idealmodel; Tr = properties.Tr, Vr = 1/properties.rhor)
    end

    ideal = _parse_ideal(ideal_data,verbose)
    #residual. it can also parse departures, that's why we pass SingleFluidResidualParam as an arg
    residual = _parse_residual(SingleFluidResidualParam,eos_data[:alphar];verbose = verbose)
    #ancillaries
    if ancillaries === nothing
        init_ancillaries = _parse_ancillaries(data[:ANCILLARIES],verbose)
        init_ancillaries.components[1] = only(_components)
    else
        init_ancillaries = init_model(ancillaries,components,ancillaries_userlocations,verbose)
    end

    references = [eos_data[:BibTeX_EOS]]

    return SingleFluid(_components,properties,init_ancillaries,ideal,residual,references)
end


"""
    SingleFluidIdeal(components;
        userlocations = String[],
        Rgas = nothing,
        verbose = false,
        coolprop_userlocations = true)

## Input parameters
- JSON data (CoolProp and teqp format)

## Input models
- `ancillaries`: a model that provides initial guesses for saturation calculations. if `nothing`, then they will be parsed from the input JSON.

## Description

Instantiates the ideal part of a single-component Empiric EoS model. `Rgas` can be used to set the value of the gas constant that is used during property calculations.

If `coolprop_userlocations` is true, then Clapeyron will try to look if the fluid is present in the CoolProp library.

The properties and ideal terms can be accessed via the `properties` and `ideal` fields respectively:

```julia-repl
julia> model = SingleFluidIdeal("water")
Ideal MultiParameter Equation of state for water:
 Lead terms: -8.3204464837497 + 6.6832105275932*τ + 3.00632*log(τ)
 Plank-Einstein terms: 5

julia> model.ideal
Ideal MultiParameter coefficients:
 Lead terms: -8.3204464837497 + 6.6832105275932*τ + 3.00632*log(τ)
 Plank-Einstein terms: 5
```

"""
function SingleFluidIdeal(components;
    userlocations = String[],
    Rgas = nothing,
    verbose = false,
    coolprop_userlocations = true,
    idealmodel = nothing,
    ideal_userlocations = String[])
    _components = format_components(components)
    single_component_check(SingleFluidIdeal,_components)
    data = get_json_data(_components;userlocations,coolprop_userlocations,verbose)    
    eos_data = first(data[:EOS])
    #properties
    properties = _parse_properties(data,Rgas,verbose)

    if idealmodel === nothing
        ideal_data = eos_data[:alpha0]
    else
        init_idealmodel = init_model(idealmodel,components,ideal_userlocations,verbose)
        ideal_data = Clapeyron.idealmodel_to_json_data(init_idealmodel; Tr = properties.Tr, Vr = 1/properties.rhor)
    end
    ideal = _parse_ideal(ideal_data,verbose)
    references = String[]
    if haskey(eos_data,:BibTeX_EOS)
        push!(references,get(eos_data,:BibTeX_EOS))
    end
    return SingleFluidIdeal(_components,properties,ideal,references)
end

function _parse_properties(data,Rgas0 = nothing, verbose = false)
    verbose && @info "Starting parsing of properties from JSON."
    #info = data[:INFO]
    eos_data_vec  = data[:EOS]
    eos_data = if eos_data_vec isa AbstractVector 
        #coolprop stores EOS field as a vector. the first one is the multiparameter #EoS
        #i did not see other examples in the CoolProp DB where they use more EoS
        first(eos_data_vec)
    else
        #this is in case we want to pass a dict directly
        eos_data_vec
    end
    st_data = data[:STATES]
    crit = st_data[:critical]
    eos_st_data = eos_data[:STATES]
    reducing = get(eos_st_data,:reducing,nothing)

    Mw = 1000*tryparse_units(get(eos_data,:molar_mass,NaN),get(eos_data,:molar_mass_units,""))
    
    T_c = tryparse_units(get(crit,:T,NaN),get(crit,:T_units,""))
    P_c = tryparse_units(get(crit,:p,NaN),get(crit,:p_units,""))
    rho_c = tryparse_units(get(crit,:rhomolar,NaN),get(crit,:rhomolar_units,""))
    if reducing !== nothing
        Tr = tryparse_units(get(reducing,:T,NaN),get(reducing,:T_units,""))
        rhor = tryparse_units(get(reducing,:rhomolar,NaN),get(reducing,:rhomolar_units,""))
    else
        Tr = T_c
        rhor = rho_c
    end
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
    a1 = 0.0 #a1
    a2 = 0.0 #a2*τ
    c0 = 0.0 #c0*log(τ)
    c1 = 0.0 #c1*τ*log(τ) (appears in one specific case)
    R0 = 0.0
    n = Float64[]
    t = Float64[]
    c = Float64[]
    d = Float64[]
    np = Float64[]
    tp = Float64[]
    n_gerg = Float64[]
    v_gerg = Float64[]
    paramtype = "ideal"
    verbose && @info "Starting parsing of $(paramtype) JSON."
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
            _a1,_a2,_c0 = _Cp0_constant_parse(id_data_i[:cp_over_R],id_data_i[:Tc],id_data_i[:T0])
            a1 += _a1
            a2 += _a2
            c0 += _c0
        elseif id_data_i[:type] == "IdealGasHelmholtzCP0PolyT"
            _T0 = id_data_i[:T0]
            _Tc = id_data_i[:Tc]
            cp_c = id_data_i[:c]
            cp_t = id_data_i[:t]
            for i in eachindex(cp_t)
                ti = cp_t[i]
                ci = cp_c[i]
                if abs(ti) <= eps(Float64) #t ≈ 0
                    #=
                    c - c * tau / tau0 + c * log(tau) ;
                    c(1 - log(tau0)) - c/tau0 * tau + c*log(tau)
                    =#
                    _a1,_a2,_c0 = _Cp0_constant_parse(ci,_Tc,_T0)
                    a1 += _a1
                    a2 += _a2
                    c0 += _c0
                elseif abs(ti + 1) <= eps(Float64) #t ≈ -1
                    #=
                    c * τ / Tc * log(τ0 / τ) + (c/Tc)*τ - (c / Tc)*τ0
                    τ*(c/Tc)*(log(τ0) - log(τ)) + (c/Tc)*τ - (c / Tc)*τ0
                    (c/Tc)*log(τ0)*τ + (c/Tc)*τ - (c / Tc)*τ0  - (c/Tc)*τ*log(τ)
                    (c/Tc)*(log(τ0) + 1)*τ - (c / Tc)*τ0  - (c/Tc)*τ*log(τ)
                    =#
                    ctc = ci/_Tc
                    tau0 = _Tc/_T0
                    a1 += -ctc*tau0
                    a2 += ctc*(log(tau0) + 1)
                    c1 += -ctc
                else
                    _a1,_a2,_npi,_tpi = _Cpi_power_parse(ci,ti,_Tc,_T0)
                    push!(np,_npi)
                    push!(tp,_tpi)
                    a1 += _a1
                    a2 += _a2
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
                _a1,_a2,_c0 = _Cp0_constant_parse(A,_Tc,_T0)
                a1 += _a1
                a2 += _a2
                c0 += _c0
            end
            n_alylee = (B,D)
            v_alylee = (C/_Tc,E/_Tc)
            append!(n_gerg,n_alylee)
            append!(v_gerg,v_alylee)
        elseif id_data_i[:type] == "IdealGasClapeyronJLGerg2008"
            append!(n_gerg,id_data_i[:n])
            append!(v_gerg,id_data_i[:v])
        elseif id_data_i[:type] == "IdealGasClapeyronJLR0"
            R0 = id_data_i[:R0]
        else
            throw(error("Ideal: $(id_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
        end
    end
    verbose && __verbose_found_json_terms(id_data)

    verbose && @info "Creating SingleFluidIdealParam from JSON."
    return SingleFluidIdealParam(a1,a2,c0,n,t,c,d,np,tp,n_gerg,v_gerg,R0)

end

function _parse_residual(out,res_data; verbose = false, Fij = 1.0)
    #polynomial y exp terms, we will separate those later
    n = Float64[]
    t = Float64[]
    d = Float64[]
    l = Float64[]
    g = Float64[]

    #gaussian terms
    n_gauss = Float64[]
    t_gauss = Float64[]
    d_gauss = Float64[]
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
    NA_C = Float64[]
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


    full = __has_extra_params(out)
    paramtype = __type_string(out)
    verbose && @info "Starting parsing of $(paramtype) JSON."

    #this is to be compatible with CoolProp departure form.
    vec_data = res_data isa AbstractVector ? res_data : (res_data,)
    for res_data_i in vec_data
        if res_data_i[:type] == "ResidualHelmholtzPower" || res_data_i[:type] == "Exponential"
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
        elseif res_data_i[:type] == "ResidualHelmholtzGaoB" && full
            append!(n_gao,res_data_i[:n])
            append!(t_gao,res_data_i[:t])
            append!(d_gao,res_data_i[:d])
            append!(eta_gao,res_data_i[:eta])
            append!(beta_gao,res_data_i[:beta])
            append!(gamma_gao,res_data_i[:gamma])
            append!(epsilon_gao,res_data_i[:epsilon])
            append!(b_gao,res_data_i[:b])
        elseif res_data_i[:type] == "ResidualHelmholtzNonAnalytic" && full
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
        elseif res_data_i[:type] == "ResidualHelmholtzAssociating"  && full
            if assoc == true
                throw(error("Residual: $(res_data_i[:type]) we only support one Associating term."))
            end
            assoc = true
            assoc_epsilonbar += res_data_i[:epsilonbar]
            assoc_kappabar += res_data_i[:kappabar]
            assoc_a += res_data_i[:a]
            assoc_m += res_data_i[:m]
            assoc_vbarn += res_data_i[:vbarn]
        elseif res_data_i[:type] == "ResidualHelmholtzGERG2008" || (res_data_i[:type] == "GERG-2008" && vec_data isa Tuple)
            #we do the conversion, as detailed in the EOS-LNG paper
            ng = res_data_i[:n]
            tg = res_data_i[:t]
            dg = res_data_i[:d]
            ηg = res_data_i[:eta]
            βg = res_data_i[:beta]
            γg = res_data_i[:gamma]
            εg = res_data_i[:epsilon]
            len = length(ηg)
            for i in 1:len
                if iszero(ηg[i]) && iszero(βg[i]) && iszero(γg[i]) && iszero(εg[i])
                    #power terms
                    push!(n,ng[i])
                    push!(t,tg[i])
                    push!(d,dg[i])
                    push!(l,0)
                    push!(g,1)
                else
                    #parse as gaussian + exponential
                    #convert to bigfloat precision, better parsing.
                    εij = big(εg[i])
                    ηij = big(ηg[i])
                    βij = big(βg[i])
                    γij = big(γg[i])
                    ω = βij*γij - ηij*εij*εij
                    if ηg[i] == 0 #simple exponential term
                        ni_new = ng[i]*exp(ω) |> Float64
                        push!(n,ni_new)
                        push!(t,tg[i])
                        push!(d,dg[i])
                        push!(l,1)
                        push!(g,βg[i])
                    else #convert to gaussian term
                        ν = 2*ηij*εij - βij
                        ξ = ν/(2*ηij)
                        ξg = ξ |> Float64
                        ni_new = ng[i]*exp(ω + ηij*ξ*ξ) |> Float64
                        push!(n_gauss,ni_new)
                        push!(t_gauss,tg[i])
                        push!(d_gauss,dg[i])
                        push!(eta,ηg[i])
                        push!(beta,0)
                        push!(gamma,0)
                        push!(epsilon,ξg)
                    end
                end
            end
        elseif res_data_i[:type] == "Gaussian+Exponential" && vec_data isa Tuple
            len = length(res_data_i[:n])
            ni = res_data_i[:n]
            ti = res_data_i[:t]
            di = res_data_i[:d]
            ηi = res_data_i[:eta]
            βi = res_data_i[:beta]
            γi = res_data_i[:gamma]
            εi = res_data_i[:epsilon]
            li = res_data_i[:l]
            for i in 1:len
                if ηi[i] == βi[i] == γi[i] == εi[i] == 0.0
                    push!(n,ni[i])
                    push!(t,ti[i])
                    push!(d,di[i])
                    push!(l,li[i])
                    push!(g,1)
                else
                    push!(n_gauss,ni[i])
                    push!(t_gauss,ti[i])
                    push!(d_gauss,di[i])
                    push!(eta,ηi[i])
                    push!(beta,βi[i])
                    push!(gamma,γi[i])
                    push!(epsilon,εi[i])
                end
            end
        else

            throw(error("$paramtype: $(res_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
        end
    end

    verbose && __verbose_found_json_terms(vec_data)

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

    verbose && @info "Creating $(string(out)) from JSON."
    if !full
        return out(Fij,_n,_t,_d,_l,_g,_η,_β,_γ,_ε)
    end

    #gao_b term
    gao_b = GaoBTerm(n_gao,t_gao,d_gao,eta_gao,beta_gao,gamma_gao,epsilon_gao,b_gao)

    #non analytical term
    na = NonAnalyticTerm(NA_A,NA_B,NA_C,NA_D,NA_a,NA_b,NA_beta,NA_n)

    #assoc terms
    assoc = Associating2BTerm(assoc_epsilonbar,assoc_kappabar,assoc_a,assoc_m,assoc_vbarn)

    #exponential term

   return SingleFluidResidualParam(_n,_t,_d,_l,_g,_η,_β,_γ,_ε;gao_b,na,assoc)
end

function __verbose_found_json_terms(data)
    res = String[]
    push!(res,"JSON types:")
    for data_i in data
        type = data_i[:type]
        additional =
        if type == "ResidualHelmholtzGERG2008" || type == "GERG-2008"
            " Converting to power, exponential and gaussian bell-shaped terms"
        elseif type == "IdealGasHelmholtzPlanckEinstein" || type == "IdealGasHelmholtzPlanckEinsteinFunctionT"
            " Converting to Generalized Plank-Einstein terms."
        elseif type == "IdealGasHelmholtzCP0Constant"
            " Converting to lead and LogTau terms."
        elseif type == "IdealGasHelmholtzCP0PolyT"
            " Converting to lead, LogTau and power terms."
        elseif type == "IdealGasHelmholtzCP0AlyLee"
            " Converting to lead, LogTau and GERG-2004 terms."
        elseif type == "Gaussian+Exponential"
            " Converting to power, exponential and gaussian bell-shaped terms."
        else
            ""
        end

    push!(res,"found $(info_color(type)) terms.$(additional)")
    end
    io = IOBuffer()
    show_pairs(io,res,quote_string = false)
    r = io |> take! |> String
    @info r
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
    rhol_anc = PolExpLiquid(_parse_ancilliary_func(rhol_data,:T_r,:reducing_value))
    return CompositeModel(["ancillaries"],gas = rhov_anc,liquid = rhol_anc,saturation = ps_anc)
end

#converting Clapeyron ideal models into SingleFluidParams

"""
    idealmodel_to_json_data(model::EoSModel;Tr = 1.0,T0 = 298.15, Vr = 1.0)

Transforms an `model::IdealModel` into a vector of dictionaries containing valid ideal multiparameter helmholtz terms.
`Tr` is the reducing temperature, `T0` is the reference temperature, `Vr` is the reducing volume.
## Example
```
julia> id = BasicIdeal(["water"])
BasicIdeal(Clapeyron.BasicIdealParam)

julia> Clapeyron.idealmodel_to_json_data(id)
1-element Vector{Dict{Symbol, Any}}:
 Dict(:T0 => 298.15, :type => "IdealGasHelmholtzCP0Constant", :cp_over_R => 2.5, :Tc => 1.0)

```
"""
function idealmodel_to_json_data(model;Tr = 1.0,T0 = 298.15,Vr = 1.0)
    if is_splittable(model)
        single_component_check(idealmodel_to_json_data,model)
    end
    return idealmodel_to_json_data(model,Tr,T0,Vr)
end

function idealmodel_to_json_data(model::BasicIdealModel,Tr,T0,Vr)
    [
        Dict(:type => "IdealGasHelmholtzLead",
            :a1 => - log(Vr) - 1.5*log(Tr) - 1,
            :a2 => 0,
        )
        Dict(
            :type => "IdealGasHelmholtzLogTau",
            :a => 1.5,
        )
    ]
end

function idealmodel_to_json_data(model::ReidIdealModel,Tr,T0,Vr)
    coeffs = model.params.coeffs[1] ./ Rgas(model)
    n = length(coeffs)
    [
        Dict(
            :type => "IdealGasHelmholtzLead",
            :a1 => - log(Vr) - log(298) + log(Tr),
            :a2 => 0,
        ),

        Dict(
            :type => "IdealGasHelmholtzLogTau",
            :a => -1,
        ),

        Dict(
            :type => "IdealGasHelmholtzCP0PolyT",
            :T0 => 298.0,
            :Tc => Tr,
            :c => [coeffs...],
            :t => collect(0:(n-1)),
        ),
    ]
end

function idealmodel_to_json_data(model::JobackIdealModel,Tr,T0,Vr)
    return idealmodel_to_json_data(ReidIdeal(model),Tr,T0,Vr)
end

function idealmodel_to_json_data(model::MonomerIdealModel,Tr,T0,Vr)
    Mwᵢ = model.params.Mw[1]*0.001
    Λᵢ = h/√(k_B*Mwᵢ/N_A) # * T^(-1/2)
    kᵢ = N_A*Λᵢ^3 #T^(-3/2)
    # monomer: a = ∑ xi * [log(xi*ki*T^-1.5/v)] - 1
    # ∑ xi * [log(xi) +  1.5*log(ki*T/v)]
    # ∑ xi * [log(xi) +  a0i(v,T)]
    #a0i(v,T) = log(ki) - log(v) + 1.5*log(Tinv)
    #a0i(v,T) = log(ki) - log(v) + log(vr) - log(vr) + 1.5*log(Tinv) + 1.5*log(Tr) - 1.5*log(Tr)
    #a0i(v,T) = log(ki) + log(vr/v) - log(vr)  - 1.5*log(Tr) + 1.5*log(Tr/Tinv)
    #a0i(v,T) = log(vr/v)  + log(ki)- log(vr) - 1.5*log(Tr) + 1.5*log(Tr/Tinv)
    #a1 = log(ki) - log(vr) - 1.5*log(Tr)
    #a2 = 1.5
    a1 = log(kᵢ) - log(Vr) - 1.5*log(Tr)
    [
        Dict(
            :type => "IdealGasHelmholtzLead",
            :a1 => a1 - 1,
            :a2 => 0.0,
            ),
        Dict(
            :type => "IdealGasHelmholtzLogTau",
            :a => 1.5,
            )
    ]
end

function idealmodel_to_json_data(model::WalkerIdealModel,Tr,T0,Vr)
    ni = model.groups.n_flattenedgroups[1]
    groups_i = model.groups.i_groups[1]
    Mwᵢ = sum(ni[k]*model.params.Mw[k] for k in groups_i)
    Nrot = model.params.Nrot.values
    Λᵢ = h/√(k_B*Mwᵢ/N_A) # * T^(-1/2)
    kᵢ = N_A*Λᵢ^3 #T^(-3/2)
    # monomer: a = ∑ xi * [log(xi*ki*T^-1.5/v)] - 1
    # ∑ xi * [log(xi) +  1.5*log(ki*T/v)]
    # ∑ xi * [log(xi) +  a0i(v,T)]
    #a0i(v,T) = log(ki) - log(v) + 1.5*log(Tinv)
    #a0i(v,T) = log(ki) - log(v) + log(vr) - log(vr) + 1.5*log(Tinv) + 1.5*log(Tr) - 1.5*log(Tr)
    #a0i(v,T) = log(ki) + log(vr/v) - log(vr)  - 1.5*log(Tr) + 1.5*log(Tr/Tinv)
    #a0i(v,T) = log(vr/v)  + log(ki)- log(vr) - 1.5*log(Tr) + 1.5*log(Tr/Tinv)
    #a1 = log(ki) - log(vr) - 1.5*log(Tr)
    #a2 = 1.5
    Nroti = sum(ni[k]*Nrot[k] for k in groups_i)/sum(ni[k] for k in groups_i)
    a1 = log(kᵢ) - log(Vr) - (1.5 + Nroti/2)*log(Tr)

    θ1 = model.params.theta1.values
    θ2 = model.params.theta2.values
    θ3 = model.params.theta3.values
    θ4 = model.params.theta4.values
    g1 = model.params.deg1.values
    g2 = model.params.deg2.values
    g3 = model.params.deg3.values
    g4 = model.params.deg4.values
    θ_vib = (θ1, θ2, θ3, θ4)
    g_vib = (g1, g2, g3, g4)

    n_pe = Float64[]
    t_pe = Float64[]
    c_pe = Float64[]
    d_pe = Float64[]
    n_power =Float64[]
    t_power = Float64[]
    for k in groups_i
        nik = ni[k]
        for v in 1:4
            gvk = g_vib[v][k]
            θvk = θ_vib[v][k]
            push!(n_power,nik*gvk*θvk/2/Tr)
            push!(t_power,1)
            push!(n_pe,gvk*nik)
            push!(t_pe,-θvk/Tr)
            push!(c_pe,1)
            push!(d_pe,-1)
        end
    end
    #res += z[i]*(
    #    sum(ni[k]*sum(g_vib[v][k]*(θ_vib[v][k]/2/T+log(1-exp(-θ_vib[v][k]/T))) for v in 1:4)
    #    for k in @groups(i)))
    [
        Dict(
            :type => "IdealGasHelmholtzLead",
            :a1 => a1 - 1,
            :a2 => 0.0,
            ),
        Dict(
            :type => "IdealGasHelmholtzLogTau",
            :a => 1.5 + Nroti/2,
            ),
        Dict(
            :type => "IdealHelmholtzPlanckEinsteinGeneralized",
            :n => n_pe,
            :t => t_pe,
            :c => c_pe,
            :d => d_pe,
            ),
        Dict(
            :type => "IdealGasHelmholtzPower",
            :n => n_power,
            :t => t_power,
            ),
    ]
end

function idealmodel_to_json_data(model::AlyLeeIdealModel,Tr,T0,Vr)
    A = model.params.A.values[1]
    B = model.params.B.values[1]
    C = model.params.C.values[1]
    D = model.params.D.values[1]
    E = model.params.E.values[1]
    F = model.params.F.values[1]
    G = model.params.G.values[1]
    H = model.params.H.values[1]
    I = model.params.I.values[1]

    [
        Dict(
            :type => "IdealGasHelmholtzCP0AlyLee",
            :Tc => Tr,
            :T0 => 298.15,
            :c => [A,B,C,D,E]
            ),
        Dict(
            :type => "IdealGasHelmholtzCP0AlyLee",
            :Tc => Tr,
            :T0 => 298.15,
            :c => [0.0,F,G,H,I]
            ),
        Dict(
        :type => "IdealGasHelmholtzLead",
        :a1 => - log(Vr),
        :a2 => 0.0,
        ),
    ]
end

function idealmodel_to_json_data(model::ShomateIdealModel,Tr,T0,Vr)
    coeffs = model.params.coeffs[1] ./ Rgas(model)
    n = length(coeffs)
    [
        Dict(
            :type => "IdealGasHelmholtzLead",
            :a1 => - log(Vr) - log(298) + log(Tr),
            :a2 => 0,
        ),

        Dict(
            :type => "IdealGasHelmholtzLogTau",
            :a => -1,
        ),

        Dict(
            :type => "IdealGasHelmholtzCP0PolyT",
            :T0 => 298.0,
            :Tc => Tr,
            :c => [coeffs...],
            :t => [0,1,2,3,-2],
        ),
    ]
end