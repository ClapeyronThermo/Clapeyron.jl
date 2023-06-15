const FIJ_TYPE = Clapeyron.PairParameter{Float64, SparseArrays.SparseMatrixCSC{Float64, Int64}}

struct EmpiricDepartureValues
    iterators::Vector{UnitRange{Int}}
    F::Float64
    n::Vector{Float64}
    t::Vector{Float64}
    d::Vector{Int}
    l::Vector{Int}
    g::Vector{Float64}
    eta::Vector{Float64}
    beta::Vector{Float64}
    gamma::Vector{Float64}
    epsilon::Vector{Float64}
    function EmpiricDepartureValues(F,n,t,d,l = Int[],g = ones(length(l)),
        eta = Float64[],beta = Float64[],gamma = Float64[], epsilon = Float64[])
        param = new(Vector{UnitRange{Int}}(undef,0),F,n,t,d,l,g,eta,beta,gamma,epsilon)
        _calc_iterators!(param)
        return param
    end
end

struct EmpiricDepartureParam <: EoSParam
    F::PairParam{Float64}
    parameters::PairParameter{EmpiricDepartureValues, SparseArrays.SparseMatrixCSC{EmpiricDepartureValues, Int64}}
end

@newmodelsimple EmpiricDeparture MultiFluidDepartureModel EmpiricDepartureParam

function EmpiricDeparture(components;userlocations = String[],verbose = false)
    params = getparams(components,["Empiric/departure/Empiric_departure_unlike.csv"],asymmetricparams = ["F","parameters"],userlocations = userlocations,verbose = verbose)
    raw_parameters = params["parameters"]
    F = params["F"]
    parsed_parameters = similar(F.values,EmpiricDepartureValues)
    s1,s2 = size(F.values)
    #parse JSON string to create EmpiricDepartureValues
    for i in 1:s1
        for j in 1:s2
            if !raw_parameters.ismissingvalues[i,j]
                Fij = F[i,j]
                if !iszero(Fij)
                    parsed_parameters[i,j] = _parse_departure(raw_parameters[i,j],Fij,verbose)
                else
                    raw_parameters.ismissingvalues[i,j] = true
                end
            end
        end
    end
    #compress
    parameters = PairParameter(raw_parameters.name,components,parsed_parameters,raw_parameters.ismissingvalues,raw_parameters.sourcecsvs,raw_parameters.sources)
    pkgparams = EmpiricDepartureParam(F,parameters)
    return EmpiricDeparture(pkgparams,verbose = verbose)
end

function multiparameter_a_res(model::EmpiricMultiFluid,V,T,z,departure::EmpiricDeparture,δ,τ,∑z = sum(z)) 
    lnδ = log(δ)
    lnτ = log(τ)
    aᵣ = multiparameter_a_res0(model,V,T,z,δ,τ,lnτ,∑z)
    _0 = zero(aᵣ)
    ℙ = departure.params.parameters.values
    ℙ_nonzeros = nonzeros(ℙ)
    isone(length(z)) && return aᵣ
    iszero(nnz(ℙ)) && return aᵣ

    @inbounds for j ∈ @comps
        zⱼ = z[j]
        for ii ∈ nzrange(ℙ, j)
            i = rows[ii]
            ℙᵢⱼ = ℙ_nonzeros[ii]
            k_pol,k_exp,k_gauss = ℙ.iterators 
            Fᵢⱼ = ℙᵢⱼ.F
            n,t,d = ℙᵢⱼ.n,ℙᵢⱼ.t,ℙᵢⱼ.d
            #strategy for storing.
            #n, t, d, gauss values, always require views
            #l, b does not require views. they are used just once.

            #Polynomial terms
            n_pol = view(n,k_pol)
            t_pol = view(t,k_pol)
            d_pol = view(d,k_pol)
            αᵣ += term_ar_pol(δ,τ,lnδ,lnτ,αᵣ,n_pol,t_pol,d_pol)

            #Exponential terms.
            if length(k_exp) != 0
                l,g = ℙᵢⱼ.l,ℙᵢⱼ.g
                n_exp = view(n,k_exp)
                t_exp = view(t,k_exp)
                d_exp = view(d,k_exp)
                αᵣ += term_ar_exp(δ,τ,lnδ,lnτ,αᵣ,n_exp,t_exp,d_exp,l,g)
            end

            #Gaussian bell-shaped terms
            
            if length(k_gauss) != 0
                η,β,γ,ε = ℙᵢⱼ.eta,ℙᵢⱼ.beta,ℙᵢⱼ.gamma,ℙᵢⱼ.epsilon
                n_gauss = view(n,k_gauss)
                t_gauss = view(t,k_gauss)
                d_gauss = view(d,k_gauss)
                αᵣ += term_ar_gauss(δ,τ,lnδ,lnτ,αᵣ,n_gauss,t_gauss,d_gauss,η,β,γ,ε)
            end
            
            aᵣ +=z[i]*zⱼ*Fᵢⱼ*aij
        end
     end
    return aᵣ
end

function _parse_departure(json_string::String,Fij::Float64,verbose = false)
    res_data = JSON3.read(json_string)
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
        elseif res_data_i[:type] == "ResidualHelmholtzExponential"
            append!(n,res_data_i[:n])
            append!(t,res_data_i[:t])
            append!(d,res_data_i[:d])
            append!(l,res_data_i[:l])
            append!(g,res_data_i[:g])
        elseif res_data_i[:type] == "ResidualHelmholtzGERG2008"
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
                    push!(n,ni_new)
                    push!(t,tg[i])
                    push!(d,dg[i])
                    push!(eta,ηg[i])
                    push!(beta,0)
                    push!(gamma,0)
                    push!(epsilon,ξg)
                end
            end
        else
            throw(error("Departure: $(res_data_i[:type]) not supported for the moment. open an issue in the repository for help."))
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
    
   return EmpiricDepartureValues(Fij,_n,_t,_d,_l,_g,_η,_β,_γ,_ε)
end