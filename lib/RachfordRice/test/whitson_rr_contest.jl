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
        push!(result,(K = ki, z = zi))
    end
    return result
end

function whitson_get_cases()
    try
        k_values = String(take!(Downloads.download("https://https://raw.githubusercontent.com/ClapeyronThermo/Rachford-Rice-Contest/503b92f1b2847c4459326841f538739bcb9d629f/k-values.csv", IOBuffer())))
        comps = String(take!(Downloads.download("https://https://raw.githubusercontent.com/ClapeyronThermo/Rachford-Rice-Contest/503b92f1b2847c4459326841f538739bcb9d629f/compositions.csv", IOBuffer())))
        return to_nt(k_values,comps)
    catch
        return nothing
    end
end

#=
function save_compatible(path::AbstractString, data; target_version=UInt8(24))
    buf = IOBuffer()
    serialize(buf, data)
    bytes = take!(buf)
    bytes[4] = target_version   # 1-based index 4 == header offset 3
    write(path, bytes)
end =#

@testset "WhitsonAS - Rachford-Rice Contest" begin
    CASES = whitson_get_cases()
    if CASES !== nothing
        test_cases(CASES)
    else
        @warn "couldn't test whitson RR-contest, no internet available."
    end
end