using FunctionMaps:
    MapStyle, IsMap, NotMap,
    functionmap, checkmap

struct MySimpleMap end
FunctionMaps.MapStyle(::Type{MySimpleMap}) = IsMap()

@testset "map interface" begin
    @test_throws ArgumentError checkmap((5,))   # a tuple is not a map
    @test MapStyle(2.0) isa IsMap
    @test MapStyle(2.0) == MapStyle(typeof(2.0))
    m = LinearMap(2)
    @test MapStyle(m) isa IsMap
    @test functionmap(m) === m
    @test checkmap(m) === m
    m2 = MySimpleMap()
    @test MapStyle(m2) isa IsMap
    @test checkmap(m2) == m2
end
