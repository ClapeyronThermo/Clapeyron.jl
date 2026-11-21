using Test, Clapeyron
struct LateParam{T} <: Clapeyron.ClapeyronParam
    value::T
end
Base.eltype(::Type{LateParam{T}}) where T = T
struct LateContainer{T} <: Clapeyron.ParametricEoSParam{T}
    p::LateParam{T}
end
@testset "late extension trait" begin
    p=LateParam(2.0)
    v=@inferred Clapeyron.build_parametric_param(LateContainer, p)
    @test v isa LateContainer{Float64}
    @test v.p===p
end
@testset "numeric promotion metadata and no mutation" begin
    for T in (Float32, Float64, BigFloat)
        names=["a", "b"]
        vals=T[1 2; 3 4]
        a=PairParam("A", names, vals)
        b=PairParam("B", names, copy(vals))
        c=PairParam("C", names, [0 0; 0 0])
        r=SingleParam("R", names, T[1, 2])
        q=SingleParam("Q", names, T[3, 4])
        p=@inferred Clapeyron.UNIFACParam(a, b, c, r, q)
        @test p isa Clapeyron.UNIFACParam{T}
        @test eltype(p.C.values)==T
        @test p.C.values==zeros(T, 2, 2)
        @test p.A.values==vals
        @test p.A.components==a.components
        @test p.A.ismissingvalues==a.ismissingvalues
        @test p.C.sources==c.sources
        @test p.C.sourcecsvs==c.sourcecsvs
        @test c.values==zeros(Int, 2, 2)
        @test eltype(c.values)==Int
    end
end
