abstract type MultiFluidModel <: EmpiricHelmholtzModel end

struct MultiFluidPropertyParam <: EoSParam
    Mw::SingleParam{Float64} #SingleParam
    rhoc::SingleParam{Float64} #SingleParam
    vc::SingleParam{Float64} #SingleParam
    Tc::SingleParam{Float64} #SingleParam
    pc::SingleParam{Float64} #SingleParam
    lb_v::SingleParam{Float64} #SingleParam
end

struct MultiFluidIdealParam <: EoSParam
    n0::SingleParam{NTuple{7, Float64}}
    theta::SingleParam{NTuple{4, Float64}}
    gamma_T::PairParam{Float64}
    gamma_v::PairParam{Float64}
    beta_T::PairParam{Float64}
    beta_v::PairParam{Float64}
end

struct MultiFluidSingleParam <: EoSParam
    n::PackedVectorSingleParam{Float64}
    t::PackedVectorSingleParam{Float64}
    d::PackedVectorSingleParam{Int}
    c::PackedVectorSingleParam{Int}
end

const FIJ_TYPE = Clapeyron.PairParameter{Float64, SparseArrays.SparseMatrixCSC{Float64, Int64}, SubArray{Float64, 1, Base.ReshapedArray{Float64, 1, SparseArrays.SparseMatrixCSC{Float64, Int64}, Tuple{Base.MultiplicativeInverses.SignedMultiplicativeInverse{Int64}}}, Tuple{StepRange{Int64, Int64}}, false}} 

struct MultiFluidPairParam <: EoSParam
    nij::PackedSparsePairParam{Float64} #SparsePairParam #done
    tij::PackedSparsePairParam{Float64} #SparsePairParam #done
    dij::PackedSparsePairParam{Int} #SparsePairParam #done
    Fij::FIJ_TYPE #SparsePairParam #done
    beta_ij::PackedSparsePairParam{Float64} #SparsePairParam
    gamma_ij::PackedSparsePairParam{Float64} #SparsePairParam
    eta_ij::PackedSparsePairParam{Float64} #SparsePairParam
    epsilon_ij::PackedSparsePairParam{Float64} #SparsePairParam
end