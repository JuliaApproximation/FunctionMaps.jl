using FunctionMaps:
    VcatMap,
    VectorProductMap,
    TupleProductMap

function test_product_map(T)
    ma = IdentityMap{T}()
    mb = interval_map(T(0), T(1), T(2), T(3))
    mc = interval_map(T(2), T(3), T(4), T(5))

    r1 = suitable_point_to_map(ma)
    r2 = suitable_point_to_map(ma)
    r3 = suitable_point_to_map(ma)
    r4 = suitable_point_to_map(ma)
    r5 = suitable_point_to_map(ma)

    m1 = productmap(ma,mb)
    test_generic_map(m1)
    @test isaffinemap(m1)
    @test !islinearmap(m1)
    @test !isconstantmap(m1)
    @test m1(SVector(r1,r2)) ≈ SVector(ma(r1),mb(r2))
    @test mapsize(m1) == (2,2)
    @test !hascanonicalmap(m1)
    @test !hasequalmap(m1)
    @test !hasequivalentmap(m1)

    m2 = productmap(m1,mb)
    test_generic_map(m2)
    @test m2(SVector(r1,r2,r3)) ≈ SVector(ma(r1),mb(r2),mb(r3))
    m3 = productmap(mb,m2)
    test_generic_map(m3)
    @test m3(SVector(r1,r2,r3,r4)) ≈ SVector(mb(r1),ma(r2),mb(r3),mb(r4))
    m4 = productmap(m1,m2)
    test_generic_map(m4)
    @test m4(SVector(r1,r2,r3,r4,r5)) ≈ SVector(m1(SVector(r1,r2))...,m2(SVector(r3,r4,r5))...)

    @test ProductMap(ma, mb) == ProductMap([ma,mb])
    @test hash(ProductMap(ma, mb)) == hash(ProductMap([ma,mb]))

    m5 = productmap(AffineMap(SMatrix{2,2,T}(1.0,2,3,4), SVector{2,T}(1,3)), LinearMap{T}(2.0))
    @test domaintype(component(m5,1)) == SVector{2,T}
    @test domaintype(component(m5,2)) == T
    test_generic_map(m5)
    x = SVector{3,T}(rand(T,3))
    @test m5(x) ≈ SVector(component(m5,1)(SVector(x[1],x[2]))...,component(m5,2)(x[3]))

    m6 = ProductMap([ma,mb])
    @test m6 isa VectorProductMap
    @test mapsize(m6) == (2,2)
    @test convert(Map{SVector{2,T}}, m6) isa VcatMap
    test_generic_map(m6)

    @test FunctionMaps.compatibleproductdims(m1, m6)
    @test composedmap(m1, m6) isa ProductMap
    @test !FunctionMaps.compatibleproductdims(m2, m5)
    @test composedmap(m2, m5) isa ComposedMap

    @test ProductMap(SA[mb,mc]) isa VcatMap
    @test ProductMap{Tuple{T,T}}(mb, mc) isa TupleProductMap

    @test components(productmap(LinearMap(1), LinearMap(one(T)))) isa Tuple{<:Map{T},<:Map{T}}

    # test a non-rectangular product map
    m_rect = AffineMap(SA[one(T) 2; 3 4; 5 6], SA[one(T),one(T),zero(T)])
    pm = productmap(m_rect, m_rect)
    @test pm isa VcatMap{T,6,4,(3,3),(2,2)}
    @test isaffinemap(pm)
    @test jacobian(pm, SA[one(T),0,0,0]) isa SMatrix{6,4,T}
    @test diffvolume(pm, SA[one(T),0,0,0]) == 24
    @test factors(pm) == components(pm)
    @test nfactors(pm) == ncomponents(pm)
    @test factor(pm, 1) == component(pm, 1)

    @test mapconstant(productmap(ConstantMap(2), ConstantMap(4))) == [2,4]

    @test VcatMap([ma,mb]) isa VcatMap
    @test VcatMap{T,2,2}([ma,mb]) isa VcatMap{T,2,2}
    @test VectorProductMap([ma,mb]) isa VectorProductMap
    @test VectorProductMap{Vector{Float64}}([ma,mb]) isa VectorProductMap{Vector{Float64}}
    @test TupleProductMap((ma,mb)) isa TupleProductMap
    @test TupleProductMap{Tuple{T,T}}((ma,mb)) isa TupleProductMap{Tuple{T,T}}
end
