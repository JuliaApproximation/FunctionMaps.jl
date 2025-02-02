function test_composite_maps()
    T = Float64
    a = T(0)
    b = T(1)
    c = T(2)
    d = T(3)
    ma = IdentityMap{T}()
    mb = interval_map(a, b, c, d)

    r = suitable_point_to_map(ma)
    m1 = ComposedMap(mb, ma)
    test_generic_map(m1)
    @test m1(r) ≈ ma(mb(r))
    m2 = ComposedMap(mb, m1)
    test_generic_map(m2)

    @test ncomponents(composedmap(m1,m1)) == 4

    m5 = ComposedMap(LinearMap(rand(T,2,2)), AffineMap(rand(T,2,2),rand(T,2)))
    test_generic_map(m5)
    @test jacobian(m5) isa ConstantMap
    @test m5[Component(1)] isa LinearMap
    @test m5[Component(2)] isa AffineMap
    @test ComposedMap(m5[Component(1:2)]...) == m5
    @test_throws BoundsError m5[Component(3)]

    m6 = multiply_map(ma,ma)
    @test m6(one(T)/2) == one(T)/4
    @test jacobian(m6) isa SumMap
    @test jacobian(m6)(one(T)) == 2
    @test jacobian(m6, one(T)) == 2

    m7 = sum_map(ma,ma)
    @test m7(one(T)) == 2
    @test jacobian(m7) isa ConstantMap
    @test jacobian(m7, one(T)) == 2
    @test sum_map() == ()
    @test sum_map(ma) == ma
    @test sum_map(ma, mb, ma) isa SumMap{T,Tuple{X,Y,Z}} where X where Y where Z
    @test sum_map(SumMap(ma, mb), ma) isa SumMap{T,Tuple{X,Y,Z}} where X where Y where Z
    @test sum_map(ma, SumMap(ma, mb)) isa SumMap{T,Tuple{X,Y,Z}} where X where Y where Z
    @test sum_map(SumMap(ma, mb), SumMap(ma,mb)) isa SumMap{T,Tuple{W,X,Y,Z}} where W where X where Y where Z
    @test FunctionMaps.SumMap{Float64}(LinearMap(2), LinearMap(2.0)) isa FunctionMaps.SumMap{Float64}
    @test allequal(map(domaintype,components(FunctionMaps.SumMap{Float64}(LinearMap(2), LinearMap(2.0)))))

    @test mapsize(ComposedMap(LinearMap(2),LinearMap(rand(T,2)),LinearMap(rand(T,2,2)))) == (2,)
    @test mapsize(ComposedMap(LinearMap(2),LinearMap(rand(T,2)))) == (2,)
    @test mapsize(ComposedMap(LinearMap(rand(T,2)),LinearMap(rand(T,2)'),LinearMap(rand(T,2)))) == (2,)
    @test mapsize(ComposedMap(LinearMap(rand(T,2,2)),LinearMap(rand(T,2)'),LinearMap(2))) == (1,2)
    @test mapsize(ComposedMap(LinearMap(one(T)),LinearMap(one(T)))) == ()

    @test composedmap() == ()
    @test composedmap(ma) == ma
    @test composedmap(ma,ma) == ma
    @test composedmap(ma,ma,ma) == ma

    @test composite_jacobian(ma) == jacobian(ma)

    @test multiply_map() == ()
    @test multiply_map(ma) == ma
    @test multiply_map(ma,ma)(2*one(T)) == 4
    @test multiply_map(ma,ma,ma)(2*one(T)) == 8
    @test FunctionMaps.mul_jacobian() == ()
    @test FunctionMaps.mul_jacobian(ma) == jacobian(ma)
    @test jacobian(multiply_map(ma, ma, ma)) isa FunctionMaps.SumMap
    
    @test ncomponents(multiply_map(multiply_map(ma,ma),multiply_map(ma,ma))) == 4
    @test multiply_map(LinearMap(2), LinearMap(2.0)) isa FunctionMaps.MulMap{Float64}
    @test FunctionMaps.MulMap{Float64}(LinearMap(2), LinearMap(2.0)) isa FunctionMaps.MulMap{Float64}
    @test allequal(map(domaintype,components(FunctionMaps.MulMap{Float64}(LinearMap(2), LinearMap(2.0)))))

    @test sum_jacobian() == ()
    @test sum_jacobian(ma) == jacobian(ma)

    @test convert_domaintype(Float64, ComposedMap(LinearMap(2), LinearMap(2))) isa Map{Float64}
    @test allequal(map(domaintype,components(convert_domaintype(Float64, ComposedMap(LinearMap(2), LinearMap(2))))))
end
