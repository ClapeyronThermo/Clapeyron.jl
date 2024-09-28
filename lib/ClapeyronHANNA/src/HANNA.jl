
struct HANNAParam <: EoSParam
    smiles::SingleParam{String}
    emb_scaled::Vector{Vector{Float64}}
    T_scaler::Function #TODO:maybe parametrize this?
    theta::Dense
    alpha::Chain
    phi::Chain
    Mw::SingleParam{Float64}
end

function Clapeyron.split_model(param::HANNAParam,splitter)
    return [Clapeyron.each_split_model(param,i) for i ∈ splitter]
end 

function Clapeyron.each_split_model(param::HANNAParam,i)
    Mw = each_split_model(param.Mw,i)
    smiles = each_split_model(param.smiles,i)
    emb_scaled = each_split_model(param.emb_scaled,i)
    return HANNAParam(smiles,emb_scaled,param.T_scaler,param.theta,param.alpha,param.phi,Mw)
end 

abstract type HANNAModel <: ActivityModel end

struct HANNA{c<:EoSModel} <: HANNAModel
    components::Array{String,1}
    params::HANNAParam
    puremodel::EoSVectorParam{c}
    references::Array{String,1}
end

export HANNA

"""
    HANNA <: ActivityModel
    HANNA(components;
    puremodel = nothing,
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `canonicalsmiles`: canonical SMILES representation of the components
- `Mw`: Single Parameter (`Float64`) (Optional) - Molecular Weight `[g/mol]`

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
Hard-Constraint Neural Network for Consistent Activity Coefficient Prediction (HANNA).
The implementation is based on 

To use the model, the package `ClapeyronHANNA` must be installed and loaded (see example below).

Recommended usage to ensure canonical smiles:
```julia
using Clapeyron, ClapeyronHANNA, ChemicalIdentifier

components = ["water","isobutanol"]
chemids = [search_chemical(c) for c in components]

model = HANNA(components,userlocations=(;Mw=[chemids[i].MW for i in 1:2],smiles=[chemids[i].smiles for i in 1:2]))
```

## References
1. Specht, T., Nagda, M., Fellenz, S., Mandt, S., Hasse, H., Jirasek, F., HANNA: Hard-Constraint Neural Network for Consistent Activity Coefficient Prediction. arXiv July 25, 2024. [10.48550/arXiv.2407.18011](https://doi.org/10.48550/arXiv.2407.18011).
"""
HANNA

default_locations(::Type{HANNA}) = ["properties/identifiers.csv", "properties/molarmass.csv"]

function HANNA(components;
        puremodel = BasicIdeal,
        userlocations = String[],
        pure_userlocations = String[],
        verbose = false)

    # Get parameters (Mw and smiles)
    params = getparams(components,default_locations(HANNA);userlocations=userlocations,ignore_headers=["dipprnumber","inchikey","cas"])

    # Load ChemBERTa model
    chembert_model, chembert_tokenizer = load_chembert()
    embs = [Vector{Float64}(undef,N_emb_chembert) for _ in eachindex(components)]
    for i in eachindex(components)
        embs[i] = chembert_model(encode(chembert_tokenizer, params["canonicalsmiles"][i])).hidden_state[:,1,1]
    end

    # Load scalers and scale embeddings
    T_scaler, emb_scaler = load_scalers()
    emb_scaled = emb_scaler.(embs)

    # Component embedding network
    w_θ = CSV.File("$DB_PATH/HANNA/theta_1_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b_θ = parse.(Float64,readlines("$DB_PATH/HANNA/theta_1_b.csv"))
    theta = Dense(N_emb_chembert => N_nodes, silu; init=(_,_)->w_θ)
    theta.bias .+= b_θ

    # Mixture embedding network
    w1_α = CSV.File("$DB_PATH/HANNA/alpha_1_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b1_α = parse.(Float64,readlines("$DB_PATH/HANNA/alpha_1_b.csv"))
    w2_α = CSV.File("$DB_PATH/HANNA/alpha_2_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b2_α = parse.(Float64,readlines("$DB_PATH/HANNA/alpha_2_b.csv"))
    alpha = Chain(
        Dense(N_nodes + 2 => N_nodes, silu; init=(_,_)->w1_α),
        Dense(N_nodes => N_nodes, silu; init=(_,_)->w2_α),
    )
    alpha[1].bias .+= b1_α
    alpha[2].bias .+= b2_α

    # Property network input
    w1_ϕ = CSV.File("$DB_PATH/HANNA/phi_1_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b1_ϕ = parse.(Float64,readlines("$DB_PATH/HANNA/phi_1_b.csv"))
    w2_ϕ = CSV.File("$DB_PATH/HANNA/phi_2_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b2_ϕ = parse.(Float64,readlines("$DB_PATH/HANNA/phi_2_b.csv"))
    phi = Chain(
        Dense(N_nodes => N_nodes, silu; init=(_,_)->w1_ϕ),
        Dense(N_nodes => 1; init=(_,_)->w2_ϕ),
    )
    phi[1].bias .+= b1_ϕ
    phi[2].bias .+= b2_ϕ
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)
    params = HANNAParam(params["canonicalsmiles"],emb_scaled,T_scaler,theta,alpha,phi,params["Mw"])
    references = String["10.48550/arXiv.2407.18011"]
    
    return HANNA(components,params,_puremodel,references)
end

function C.excess_gibbs_free_energy(model::HANNAModel,p,T,z)
    x = z ./ sum(z)

    # Scale input (T and embs)
    T_s = model.params.T_scaler(T)
    # Fine tuning of the component embeddings
    θ_i = model.params.theta.(model.params.emb_scaled)
    gE = zero(Base.promote_eltype(T_s,x))
    n = length(model)
    for i in 1:n 
        for j in (i+1):n
            # Calculate cosine similarity and distance between the two components
            cosine_sim_ij = cosine_similarity(θ_i[i],θ_i[j])
            cosine_dist_ij = 1.0 - cosine_sim_ij

            # Concatenate embeddings with T and x
            c_i = vcat.(T_s,[x[i],x[j]],[θ_i[i],θ_i[j]])
            α_i = model.params.alpha.(c_i)
            c_mix = sum(α_i)
            gE_NN = model.params.phi(c_mix)[1]

            # Apply cosine similarity adjustment
            gE += x[i]*x[j]*gE_NN*cosine_dist_ij
        end
    end
    return gE * Rgas(model) * T * sum(z)
end
