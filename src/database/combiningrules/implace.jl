### kij_mix

function kij_mix!(f::F,out::PairParameter,::Nothing) where F
    N = length(out.components)
    k = FillArrays.Zeros(N,N)
    out_missing = out.ismissingvalues
    kij_mix!(f,out.values,k,out_missing)
    #if kij is missing, then the output values should be the same as the input values.
    #no missing prop has to be done
    return out
end

function kij_mix!(f::F,out::PairParameter) where F
    return kij_mix!(f,out,nothing)
end

function kij_mix!(f::F,out::PairParameter,K::PairParameter) where F
    out_missing = out.ismissingvalues
    kij_mix!(f,out.values,K.values,out_missing)
    #should consider the two.
    out_missing .= out_missing .& K.ismissingvalues
    #but diagonals are all non-missing, by default:
    diagvalues(out_missing) .= false
    return out
end

function kij_mix!(f::F,out::PairParameter,K::AbstractMatrix) where F
    out_missing = out.ismissingvalues
    kij_mix!(f,out.values,K.values,out_missing)
    #should consider the two.
    out_missing .= out_missing .& K.ismissingvalues
    #but diagonals are all non-missing, by default:
    diagvalues(out_missing) .= false
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

function pair_mix!(f::F,out::PairParameter,Q::PairParameter) where F
    out_missing = out.ismissingvalues
    pair_mix!(f,out.values,Q.values,out_missing)
    #consider the two here:
    out_missing .= out_missing .& Q.ismissingvalues
    #but diagonals are all non-missing, by default:
     for i in diagind(out_missing)
        out_missing[i] = false
    end
    return out
end

function pair_mix!(f::F,out::PairParameter,Q::SingleParameter) where F
    out_missing = out.ismissingvalues
    pair_mix!(f,out.values,Q.values,out_missing)
    #Q is single, so all diagonal values are zero, not set
    #out_missing .= out_missing .& Q.ismissingvalues
    #but diagonals are all non-missing, by default:
    diagvalues(out_missing) .= false
    return out
end


#f: function of the form f(pi,pk,k)
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
                q_ji =  q[j,i]
                p_ji = f(p_j,p_i,q_j,q_i,q_ji)
                p[j,i] = p_ji
            end
        end
    end
    return p
end