include("HANNA_utils.jl")

struct HANNAParam <: EoSParam
    emb_scaled::Vector{Vector{Float64}}
    T_scaler::Function
    theta::Dense
    alpha::Chain
    phi::Chain
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
    puremodel = BasicIdeal(),
    userlocations = String[],
    pure_userlocations = String[],
    verbose = false)

## Input parameters
- `smiles`: SMILES representation of the components
- `Mw`: Single Parameter (`Float64`) (Optional) - Molecular Weight `[g/mol]`

## Input models
- `puremodel`: model to calculate pure pressure-dependent properties

## Description
Hard-Constraint Neural Network for Consistent Activity Coefficient Prediction (HANNA)

Recommended usage to ensure using canonical smiles:
```julia
using Clapeyron, ChemicalIdentifier

components = ["water","isobutanol"]
chemids = [search_chemical(c) for c in components]

model = HANNA(components,userlocations=(;Mw=[chemids[i].MW for i in 1:2],smiles=[chemids[i].smiles for i in 1:2]))
```

## References
1 Specht, T., Nagda, M., Fellenz, S., Mandt, S., Hasse, H., Jirasek, F., HANNA: Hard-Constraint Neural Network for Consistent Activity Coefficient Prediction. arXiv July 25, 2024. [10.48550/arXiv.2407.18011](https://doi.org/10.48550/arXiv.2407.18011).

"""
HANNA

default_locations(::Type{HANNA}) = ["properties/identifiers.csv", "properties/molarmass.csv"]

function HANNA(components;
        puremodel = BasicIdeal(),
        userlocations = String[],
        pure_userlocations = String[],
        verbose = false)

    # Get parameters (Mw and smiles)
    params = getparams(components,default_locations(HANNA);userlocations=userlocations,ignore_headers=["dipprnumber","inchikey","cas"])

    # Load ChemBERTa model
    chembert_model, chembert_tokenizer = load_chembert()
    embs = [Vector{Float64}(undef,N_emb_chembert) for _ in eachindex(components)]
    for i in eachindex(components)
        embs[i] = chembert_model(encode(chembert_tokenizer, params["smiles"][i])).hidden_state[:,1,1]
    end

    # Load scalers and scale embeddings
    T_scaler, emb_scaler = load_scalers()
    emb_scaled = emb_scaler.(embs)

    # Load HANNA model parameters
    DB = SHORT_PATHS["DB"]

    # Component embedding network
    w_θ = CSV.File("$DB/Activity/HANNA/theta_1_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b_θ = parse.(Float64,readlines("$DB/Activity/HANNA/theta_1_b.csv"))
    theta = Dense(N_emb_chembert => N_nodes, silu; init=(_,_)->w_θ)
    theta.bias .+= b_θ

    # Mixture embedding network
    w1_α = CSV.File("$DB/Activity/HANNA/alpha_1_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b1_α = parse.(Float64,readlines("$DB/Activity/HANNA/alpha_1_b.csv"))
    w2_α = CSV.File("$DB/Activity/HANNA/alpha_2_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b2_α = parse.(Float64,readlines("$DB/Activity/HANNA/alpha_2_b.csv"))
    alpha = Chain(
        Dense(N_nodes + 2 => N_nodes, silu; init=(_,_)->w1_α),
        Dense(N_nodes => N_nodes, silu; init=(_,_)->w2_α),
    )
    alpha[1].bias .+= b1_α
    alpha[2].bias .+= b2_α

    # Property network input
    w1_ϕ = CSV.File("$DB/Activity/HANNA/phi_1_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b1_ϕ = parse.(Float64,readlines("$DB/Activity/HANNA/phi_1_b.csv"))
    w2_ϕ = CSV.File("$DB/Activity/HANNA/phi_2_w.csv";header=0,types=Float64,delim=',') |> CSV.Tables.matrix
    b2_ϕ = parse.(Float64,readlines("$DB/Activity/HANNA/phi_2_b.csv"))
    phi = Chain(
        Dense(N_nodes => N_nodes, silu; init=(_,_)->w1_ϕ),
        Dense(N_nodes => 1; init=(_,_)->w2_ϕ),
    )
    phi[1].bias .+= b1_ϕ
    phi[2].bias .+= b2_ϕ
    
    _puremodel = init_puremodel(puremodel,components,pure_userlocations,verbose)

    params = HANNAParam(emb_scaled,T_scaler,theta,alpha,phi)
    references = String["10.48550/arXiv.2407.18011"]
    
    return HANNA(components,params,_puremodel,references)
end

function excess_gibbs_free_energy(model::HANNAModel,p,T,x)
    # x = z./sum(z)
    x ./ sum(x)

    # Scale input (T and embs)
    T_s = model.params.T_scaler(T)
    
    # Fine tuning of the component embeddings
    θ_i = model.params.theta.(model.params.emb_scaled)

    # Calculate cosine similarity and distance between the two components
    cosine_sim = cosine_similarity(θ_i[1],θ_i[2])
    cosine_dist = 1.0-cosine_sim

    # Concatenate embeddings with T and x
    c_i = vcat.(T_s,x,θ_i)
    α_i = model.params.alpha.(c_i)
    c_mix = sum(α_i)
    gE_NN = model.params.phi(c_mix)[1]

    # Apply cosine similarity adjustment
    return gE_NN * prod(x) * cosine_dist * Rgas(model) * T
end


function activity_coefficient(model::HANNAModel,p,T,z)
    gE = excess_gibbs_free_energy(model,p,T,z) / Rgas(model) / T
    dgEdx1 = ForwardDiff.derivative(x -> excess_gibbs_free_energy(model,p,T,[x, 1-x]) / Rgas(model) / T, z[1])
    lnγ1x = gE .+ z[2].*dgEdx1
    lnγ2x = gE .- z[1].*dgEdx1
    return exp.([lnγ1x, lnγ2x])
end