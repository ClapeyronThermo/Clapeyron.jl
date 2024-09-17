using Transformers.HuggingFace, Transformers.TextEncoders, Flux

# Constants
const N_emb_chembert = 384
const N_nodes = 96

"""
    load_chembert()

Load ChemBERTa model from HuggingFace (`DeepChem/ChemBERTa-77M-MTR`).
"""
function load_chembert(;name="DeepChem/ChemBERTa-77M-MTR", max_length=128)
    config = load_config(name)
    myconfig = HuggingFace.HGFConfig(config; max_length=max_length)
    return load_model(name; config=myconfig), load_tokenizer(name; config=myconfig)
end

function load_scalers()
    data = CSV.File("$(SHORT_PATHS["DB"])/Activity/HANNA/scaler_data.csv";header=1,types=Float64,delim=',') |> CSV.Tables.matrix
    (u,s,k) = [data[:,i][:] for i in 1:3]

    T_scaler(x::Float64) = (x - u[1]) / s[1] *k[1]
    emb_scaler(x) = (x .- u[2:end]) ./ s[2:end] .* k[2:end]

    return T_scaler, emb_scaler
end

silu(x) = @. x/(1+exp(-x))

# Cosine similarity
function cosine_similarity(x1,x2;eps=1e-8)
    res = sum(x1.*x2)./(max.(sqrt.(sum(x1.^2)),eps).*max.(sqrt.(sum(x2.^2)),eps))
end
