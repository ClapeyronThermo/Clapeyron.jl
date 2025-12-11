## utils
SingleOrPair_values(x) = x
SingleOrPair_values(x::SingleOrPair) = x.values

### kij_mix

"""
    kij_mix!(f,out)

Inplace version of [`kij_mix`](@ref)
"""
function kij_mix!(f::F,out) where F
    return kij_mix!(f,out,nothing)
end

function kij_mix!(f::F,out,::Nothing) where F
    _out = SingleOrPair_values(out)
    N = LinearAlgebra.checksquare(_out)
    k = FillArrays.Zeros(N,N)

    is_param = out isa SingleOrPair

    if is_param
        out_missing = out.ismissingvalues
        kij_mix!(f,_out,k,out_missing)
    else
        kij_mix!(f,_out,k,FillArrays.Fill(true,(N,N)))
    end
    #if kij is missing, then the output values should be the same as the input values.
    #no missing prop has to be done
    return out
end


function kij_mix!(f::F,out,K) where F
    is_param = out isa SingleOrPair
    _out = SingleOrPair_values(out)
    N = LinearAlgebra.checksquare(_out)
    _K = SingleOrPair_values(K)

    if !is_param
        kij_mix!(f,_out,_K,FillArrays.Fill(true,(N,N)))
    else
        out_missing = out.ismissingvalues  
        if K isa PairParameter
            kij_mix!(f,_out,_K,out_missing)
            #missing propagation should consider the two inputs.
            out_missing .= out_missing .& K.ismissingvalues
            #but diagonals are all non-missing, by default:
            diagvalues(out_missing) .= false
        else
            #when passing a K::AbstractMatrix, we assume that all values of k are specified.
            #that means we want to calculate all kij interactions
            out_missing .= true
            kij_mix!(f,out.values,K,out_missing)
            out_missing .= false
        end
    end
    return out
end

#f: function of the form f(pi,pk,k)
#p: property matrix, symmetric with set diagonal values and maybe-set symmetric non-diagonal entries
#K: kij matrix, symmetric, with zero diagonal
#B: bool matrix, indicates which of the non diagonal entries of p are not set.
function kij_mix!(f::F,p::AbstractMatrix,K::AbstractMatrix,B::AbstractMatrix) where F
    N = LinearAlgebra.checksquare(p)
    for j ∈ 1:N
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
    _out = SingleOrPair_values(out)
    q = SingleOrPair_values(Q)
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