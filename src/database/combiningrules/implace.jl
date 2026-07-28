#utils
function kij_to_mat(out,::Nothing)
    N = size(raw_values(out),1)
    K = FillArrays.Zeros(N,N)
   return K
end

kij_to_mat(out,K) = raw_values(K)

function bool_to_mat(out,::Nothing)
    N = size(raw_values(out),1)
    B = FillArrays.Fill(true,N,N)
   return B
end

bool_to_mat(out,B) = raw_values(B)



### kij_mix

"""
    kij_mix!(f,out)

Implace version of [`kij_mix`](@ref)
"""
kij_mix!(f::F,out) where F = kij_mix!(f,out,nothing)

kij_mix!(f::F,out,K) where F = _kij_mix!(f,raw_values(out),bool_to_mat(out,nothing),kij_to_mat(out,K))

function kij_mix!(f::F,out::PairParameter,::Nothing) where F
    _kij_mix!(f,raw_values(out),out.ismissingvalues,kij_to_mat(out,nothing))
    #if kij is missing, then the output values should be the same as the input values.
    #no missing propagation has to be done
    return out
end

function kij_mix!(f::F,out::PairParameter,K::PairParameter) where F
    out_missing = out.ismissingvalues
    _kij_mix!(f,raw_values(out),out.ismissingvalues,raw_values(K))
    #missing propagation should consider the two inputs.
    out_missing .= out_missing .& K.ismissingvalues
    #but diagonals are all non-missing, by default:
    diagvalues(out_missing) .= false
    return out
end

function kij_mix!(f::F,out::PairParameter,K) where F
    #when passing a K::AbstractMatrix, we assume that all values of k are specified.
    #that means we want to calculate all kij interactions
    out_missing = out.ismissingvalues
    out_missing .= true
    _kij_mix!(f,raw_values(out),out.ismissingvalues,kij_to_mat(out,K))
    out_missing .= false
    return out
end

kij_mix!(f,p,B,K) = _kij_mix!(f,p,B,K)

#f: function of the form f(pi,pk,k)
#p: property matrix, symmetric with set diagonal values and maybe-set symmetric non-diagonal entries
#K: kij matrix, symmetric, with zero diagonal
#B: bool matrix, indicates which of the non diagonal entries of p are not set.
@inline function _kij_mix!(f::F,p::AbstractMatrix{T1},B::AbstractMatrix{Bool},K::AbstractMatrix{T2}) where {F,T1,T2}
    @boundscheck begin
        n1 = LinearAlgebra.checksquare(K)
        n2 = LinearAlgebra.checksquare(B)
        n3 = LinearAlgebra.checksquare(p)
        @assert n1 == n2 == n3
    end

    N = LinearAlgebra.size(p,1)
    @inbounds for j ∈ 1:N
        p_j = p[j,j]
        for i ∈ 1:N
            if B[j,i]
                p_i = p[i,i]
                p_ji = f(p_i,p_j,K[j,i])
                p[j,i] = p_ji
            end
        end
    end
    return p
end

## pair_mix!
"""
    pair_mix!(f,out,Q)

Inplace version of [`pair_mix`](@ref)
"""
function pair_mix!(f::F,out,Q) where F
    _out = raw_values(out)
    q = raw_values(Q)
    is_param = out isa SingleOrPair

    #fail if the mixing method requires using qij, but it is not provided.
    q isa AbstractVector && __requires_qij(f) && __qij_error(f)

    if !is_param
        pair_mix!(f,_out,q,FillArrays.Fill(true,size(_out)))
    else
        out_missing = out.ismissingvalues
        pair_mix!(f,_out,q,out_missing)
        if ndims(q) == 2
            out_missing .= out_missing .& Q.ismissingvalues
        end
        #but diagonals are all non-missing, by default:
        for i in diagind(out_missing)
            out_missing[i] = false
        end
    end
    return out
end

#dispatch for single vectors, so q[i] -> q[i,i]

#f: function of the form f(pi,pk,qi,qj,k)
#p: property matrix, symmetric with set diagonal values and maybe-set symmetric non-diagonal entries
#Q: qij matrix, symmetric
#B: bool matrix, indicates which of the non diagonal entries of p are not set.
function pair_mix!(f,p::AbstractMatrix,q::AbstractMatrix,B::AbstractMatrix)
    N = LinearAlgebra.checksquare(p)
    for j ∈ 1:N
        p_j = p[j,j]
        q_j = q[j,j]
        for i ∈ 1:N
            if B[j,i]
                p_i = p[i,i]
                q_i = q[i,i]
                q_ji = q[j,i]
                p_ji = f(p_j,p_i,q_j,q_i,q_ji)
                p[j,i] = p_ji
            end
        end
    end
    return p
end

function mirror_pair!(p::AbstractMatrix,missing_matrix::AbstractMatrix{Bool},f)
    s1,s2 = size(p)
    for i in 1:s2
        for j in 1:s1
            if !missing_matrix[i,j]
                if missing_matrix[j,i]
                    pji = f(p[i,j])
                    p[j,i] = pji
                    missing_matrix[j,i] = false
                end
            end
        end
    end
end