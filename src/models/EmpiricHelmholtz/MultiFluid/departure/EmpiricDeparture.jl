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

#for showing in SparseMatrix context.
Base.zero(x::EmpiricDepartureValues) = zero(typeof(x))
Base.zero(::Type{EmpiricDepartureValues}) = EmpiricDepartureValues(0.,Float64[],Float64[],Int[],Int[])

function Base.show(io::IO,x::EmpiricDepartureValues)
    print(io,"aij(")
    k_pol,k_exp,k_gauss = x.iterators
    l_pol,l_exp,l_gauss = length(k_pol),length(k_exp),length(k_gauss)
    l_pol != 0 && print(io,"pol=$l_pol")
    l_exp != 0 && print(io,"exp=$l_exp")
    l_gauss != 0 && print(io,"gauss=$l_gauss")
    print(io,")")
end

struct EmpiricDepartureParam <: EoSParam
    F::PairParam{Float64}
    parameters::PairParameter{EmpiricDepartureValues, SparseArrays.SparseMatrixCSC{EmpiricDepartureValues, Int64}}
end

@newmodelsimple EmpiricDeparture MultiFluidDepartureModel EmpiricDepartureParam

"""
GEDeparture <: MultiFluidDepartureModel
    GEDeparture(components;
    activity = UNIFAC,
    userlocations=String[],
    verbose=false)

## Input parameters
none
- `F`: Pair Parameter (`Float64`) - binary interaction parameter (no units)
- `parameters`: Pair Parameter (`String`) - JSON data containing the departure terms for the binary pair

## Description

Departure that uses empiric departure functions:

```
aᵣ = ∑xᵢaᵣᵢ(δ,τ) + Δa
Δa = ∑xᵢxⱼFᵢⱼaᵣᵢⱼ(δ,τ)

aᵣᵢⱼ = ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)*τ^(tᵢⱼ₋ₖ) + 
    ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)τ^(tᵢⱼ₋ₖ)*exp(-gᵢⱼ₋ₖδ^lᵢⱼ₋ₖ) +
    ∑nᵢⱼ₋ₖδ^(dᵢⱼ₋ₖ)τ^(tᵢⱼ₋ₖ)*exp(ηᵢⱼ₋ₖ(δ-εᵢⱼ₋ₖ)^2 + βᵢⱼ₋ₖ(τ-γᵢⱼ₋ₖ)^2)

```


"""
function EmpiricDeparture(components;userlocations = String[],verbose = false)
    params = getparams(components,["Empiric/departure/empiric_departure_unlike.csv"],asymmetricparams = ["F","parameters"],userlocations = userlocations,verbose = verbose)
    raw_parameters = params["parameters"]
    F = params["F"]
    s1,s2 = size(F.values)
    𝕊 = sparse((!).(raw_parameters.ismissingvalues))
    dropzeros!(𝕊)
    parsed_parameters = SparseMatrixCSC{EmpiricDepartureValues,Int}(𝕊.m, 𝕊.n, 𝕊.colptr, 𝕊.rowval, similar(𝕊.nzval,EmpiricDepartureValues))
    #parse JSON string to create EmpiricDepartureValues
    for i in 1:s1
        for j in 1:s2
            if !raw_parameters.ismissingvalues[i,j]
                Fij = F[i,j]
                if !iszero(Fij)
                    parsed_parameters[i,j] = _parse_residual(EmpiricDepartureValues,raw_parameters[i,j];verbose,Fij)
                else
                    #raw_parameters.ismissingvalues[i,j] = true
                end
            end
        end
    end
    #compress
    parameters = PairParameter(raw_parameters.name,components,parsed_parameters,raw_parameters.ismissingvalues,raw_parameters.sourcecsvs,raw_parameters.sources)    
    pkgparams = EmpiricDepartureParam(F,parameters)
    return EmpiricDeparture(pkgparams,verbose = verbose)
end

function multiparameter_a_res(model::MultiFluid,V,T,z,departure::EmpiricDeparture,δ,τ,∑z = sum(z)) 
    lnδ = log(δ)
    lnτ = log(τ)
    aᵣ = multiparameter_a_res0(model,V,T,z,δ,τ,lnδ,lnτ,∑z)
    _0 = zero(aᵣ)
    ℙ = departure.params.parameters.values
    ℙ_nonzeros = nonzeros(ℙ)
    rows = rowvals(ℙ)
    isone(length(z)) && return aᵣ
    iszero(nnz(ℙ)) && return aᵣ
    Δa = zero(aᵣ)
    @inbounds for j ∈ @comps
        zⱼ = z[j]
        for ii ∈ nzrange(ℙ, j)
            i = rows[ii]
            ℙᵢⱼ = ℙ_nonzeros[ii]
            Δaᵢⱼ = zero(Δa)
            k_pol,k_exp,k_gauss = ℙᵢⱼ.iterators 
            Fᵢⱼ = ℙᵢⱼ.F
            n,t,d = ℙᵢⱼ.n,ℙᵢⱼ.t,ℙᵢⱼ.d
            #strategy for storing.
            #n, t, d, gauss values, always require views
            #l, b does not require views. they are used just once.

            #Polynomial terms
            n_pol = view(n,k_pol)
            t_pol = view(t,k_pol)
            d_pol = view(d,k_pol)
            Δaᵢⱼ += term_ar_pol(δ,τ,lnδ,lnτ,Δaᵢⱼ,n_pol,t_pol,d_pol)
            #Exponential terms.
            if length(k_exp) != 0
                l,g = ℙᵢⱼ.l,ℙᵢⱼ.g
                n_exp = view(n,k_exp)
                t_exp = view(t,k_exp)
                d_exp = view(d,k_exp)
                Δaᵢⱼ += term_ar_exp(δ,τ,lnδ,lnτ,Δaᵢⱼ,n_exp,t_exp,d_exp,l,g)
            end

            #Gaussian bell-shaped terms
            if length(k_gauss) != 0
                η,β,γ,ε = ℙᵢⱼ.eta,ℙᵢⱼ.beta,ℙᵢⱼ.gamma,ℙᵢⱼ.epsilon
                n_gauss = view(n,k_gauss)
                t_gauss = view(t,k_gauss)
                d_gauss = view(d,k_gauss)
                Δaᵢⱼ += term_ar_gauss(δ,τ,lnδ,lnτ,Δaᵢⱼ,n_gauss,t_gauss,d_gauss,η,β,γ,ε)
            end
            Δa +=z[i]*zⱼ*Fᵢⱼ*Δaᵢⱼ
        end
     end
    return aᵣ + Δa/(∑z*∑z)
end

export EmpiricDeparture