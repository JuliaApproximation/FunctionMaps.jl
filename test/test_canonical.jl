using FunctionMaps:
    hascanonicalmap,
    hasequalmap,
    equivalentmap,
    hasequivalentmap,
    tofunctionmap,
    equalmap

struct MyCanonicalType <: FunctionMaps.CanonicalType end

function test_canonical()
    m = LinearMap(2.0)
    @test canonicalmap(m) == m
    @test !hascanonicalmap(m)
    @test canonicalmap(MyCanonicalType(), m) == m
    @test !hascanonicalmap(MyCanonicalType(), m)

    @test tofunctionmap(m) == m
    @test tofunctionmap(MapRef(m)) == m

    @test canonicalmap(FunctionMaps.Equal(), 5.0) isa Map
    @test hasequalmap(5.0)
    @test equalmap(5.0) isa ScalarLinearMap{Float64}
    @test isequalmap(LinearMap(5), 5.0)

    m1 = LinearMap(2.0)
    m2tuple = FunctionMaps.TupleProductMap(m1, m1)
    @test hasequivalentmap(m2tuple)
    @test equivalentmap(m2tuple) isa FunctionMaps.VcatMap
end
