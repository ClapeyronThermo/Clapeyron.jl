function get_k_mean(p::PairParam)
    return get_k_mean(p.values)
end

function get_k_geomean(p::PairParam)
    return get_k_geomean(p.values)
end

function get_k_powmean(p::PairParam,n = 2)
    return get_k_powmean(p.values, n)
end

function get_k_mean3(p::PairParam)
    return get_k_mean3(p.values)
end

function get_k_mean(p::AbstractMatrix)
    k = similar(p)
    k .= 0
    n = LinearAlgebra.checksquare(p)
    for i in 1:n
        p_i = p[i,i]
        for j in 1:n
            p_j = p[j,j]
            #kij = 0.5*(pi + pj)*(1-kij)
            k[i,j] = 1 - 2*p[i,j]/(p_i + p_j)
        end
    end
    return k
end

function get_k_geomean(p::AbstractMatrix)
    k = similar(p)
    k .= 0
    n = LinearAlgebra.checksquare(p)
    for i in 1:n
        p_i = p[i,i]
        for j in 1:n
            p_j = p[j,j]
            #pij = sqrt(pi * pj)*(1-kij)
            k[i,j] = 1 - p[i,j]/sqrt(p_i*p_j)
        end
    end
    return k
end

function get_k_powmean(p::AbstractMatrix,n=2)
    k = similar(p)
    k .= 0
    nn = LinearAlgebra.checksquare(p)
    for i in 1:nn
        p_i = p[i,i]
        for j in 1:nn
            p_j = p[j,j]
            #pij = (1-k)*(0.5*(p_i^n + p_j^n))^(1/n)
            k[i,j] = 1 - p[i,j]/((0.5*(p_i^n + p_j^n))^(1/n))
        end
    end
    return k
end

function get_k_mean3(p::AbstractMatrix)
    k = similar(p)
    k .= 0
    nn = LinearAlgebra.checksquare(p)
    for i in 1:nn
        p_i = p[i,i]
        for j in 1:nn
            p_j = p[j,j]
            #pij = (1-k)*(0.5*(cbrt(p_i) + cbrt(p_j)))^3
            k[i,j] = 1 - p[i,j]/((0.5*(cbrt(p_i) + cbrt(p_j)))^3)
        end
    end
    return k
end