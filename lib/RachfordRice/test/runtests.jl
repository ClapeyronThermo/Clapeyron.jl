#=
function to_nt(K,comps)
    result_type = @NamedTuple{nc::Int,K::Vector{Float64},z::Vector{Float64}}
    result = result_type[]

    cache = Vector{Float64}(undef,105)
    nc = comps[:Nc]
    for i in eachindex(nc)
        ki_rows = K[i]
        j = 0
        nci = nc[i]
        for Kij in ki_rows
            j += 1
            cache[j] = Kij
            j == nci && break
        end
        ki = copy(@view(cache[1:nci]))
        j = 0
        zi_rows = comps[i]
        for zij in zi_rows
            j += 1
            cache[j] = zij
            j == (nci+1) && break
        end

        zi = copy(@view(cache[2:nci+1]))
        push!(result,(nc = nci,K = ki, z = zi))
    end
    return result
end

function save_compatible(path::AbstractString, data; target_version=UInt8(24))
    buf = IOBuffer()
    serialize(buf, data)
    bytes = take!(buf)
    bytes[4] = target_version   # 1-based index 4 == header offset 3
    write(path, bytes)
end =#


@testset "WhitsonAS - Rachford-Rice Contest" begin

    failures = []
    for case in CASES
        K,z = case.K,case.z
        RR = RachfordRiceResidual(K,z)
        result = RachfordRice.rachfordrice(RR,K_tol = 0.0)
        L,V = result.L, result.V
        x = liquid_compositions(RR,result)
        y = vapour_compositions(RR,result)
        all_converged,severities, convergences = RachfordRice.RachfordRiceConvergence.check_convergence(y, x, V, L, z, K)
        if !all_converged
            push!(failures,(RR,result,severities, convergences))
        end
    end
end
