using FunctionMaps: ScalarAffineMap,
    VectorAffineMap,
    StaticAffineMap,
    GenericAffineMap,
    ScalarLinearMap,
    VectorLinearMap,
    StaticLinearMap,
    GenericLinearMap,
    ScalarTranslation,
    VectorTranslation,
    StaticTranslation,
    GenericTranslation,
    ProductMap,
    TupleProductMap, VcatMap, VectorProductMap,
    WrappedMap,
    interval_map, multiply_map,
    SumMap, sum_map,
    composedmap, ComposedMap,
    composite_jacobian, sum_jacobian,
    CartToPolarMap, PolarToCartMap,
    UnitCircleMap, AngleMap, UnitDiskMap,
    VectorToComplex


maps_to_test(T) = [
    IdentityMap{T}(),
    IdentityMap{Vector{T}}(10),
    ConstantMap{T}(one(T)),
    ConstantMap{T}(SVector{2,T}(1,2)),
    ZeroMap{T}(),
    UnityMap{T}(),
    AffineMap(T(1.2), T(2.4)), # scalar map
    AffineMap(-T(1.2), T(2.4)), # scalar map with negative A
    AffineMap(randvec(T, 2, 2), randvec(T, 2)), # static map
    AffineMap(randvec(T, 3, 2), randvec(T, 3)), # static map, rectangular
    AffineMap(rand(T, 2, 2), rand(T, 2)), # vector map
    AffineMap(rand(T, 3, 2), rand(T, 3)), # vector map, rectangular
    AffineMap(LinearAlgebra.I, one(T)/2),  # use UniformScaling object as A
    AffineMap(LinearAlgebra.I, randvec(T,2)),  # use UniformScaling object as A
    GenericAffineMap(randvec(T, 2, 2), randvec(T, 2)),
    GenericAffineMap(T(1.2), randvec(T, 2)),
    GenericAffineMap(randvec(T, 3, 2), randvec(T, 3)),
    Translation(randvec(T, 3)),
    LinearMap(randvec(T, 2, 2)),
    LinearMap(randvec(T, 2)),
    AffineMap(5.0, 2.0) ∘ VectorToComplex{T}() ∘ UnitCircleMap{T}(),
    LinearMap(SMatrix{2,2}(1,2,3,T(4))) ∘ CartToPolarMap{T}() ∘ LinearMap(SMatrix{2,2}(1,2,3,T(4))),
    # Interval{Any}(0.0, 1.0)
]

randvec(T,n) = SVector{n,T}(rand(n))
randvec(T,m,n) = SMatrix{m,n,T}(rand(m,n))

suitable_point_to_map(m::Map) = suitable_point_to_map(m, domaintype(m))

suitable_point_to_map(m::Map, ::Type{SVector{N,T}}) where {N,T} = SVector{N,T}(rand(N))
suitable_point_to_map(m::Map, ::Type{T}) where {T<:Number} = rand(T)
suitable_point_to_map(m::Map, ::Type{<:AbstractVector{T}}) where {T} = rand(T, mapsize(m,2))

suitable_point_to_map(m::ProductMap) =
    map(suitable_point_to_map, components(m))
suitable_point_to_map(m::VcatMap{T,M,N}) where {T,M,N} =
    SVector{N,T}(rand(T,N))

suitable_point_to_map(::CartToPolarMap{T}) where {T} = randvec(T,2)
suitable_point_to_map(::PolarToCartMap{T}) where {T} = randvec(T,2)

widertype(T) = widen(T)
widertype(::Type{SVector{N,T}}) where {N,T} = SVector{N,widen(T)}
widertype(::Type{Vector{T}}) where {T} = Vector{widen(T)}
widertype(::Type{Tuple{A}}) where {A} = Tuple{widen(A)}
widertype(::Type{Tuple{A,B}}) where {A,B} = Tuple{widen(A),widen(B)}
widertype(::Type{Tuple{A,B,C}}) where {A,B,C} = Tuple{widen(A),widen(B),widen(C)}
widertype(::Type{NTuple{N,T}}) where {N,T} = NTuple{N,widen(T)}

issquarematrix(A) = false
issquarematrix(A::AbstractArray) = size(A,1)==size(A,2)


function test_maps()
    @testset "generic functionality" begin
        test_generic_functionality()
        test_canonical()
    end

    @testset "generic map tests" begin
        generic_map_tests(Float64)
        generic_map_tests(BigFloat)
    end

    # Test special maps
    @testset "identity map" begin
        test_identity_map(Float64)
        test_identity_map(BigFloat)
    end
    @testset "basic maps" begin
        test_basic_maps(Float64)
        test_basic_maps(BigFloat)
    end
    @testset "affine maps" begin
        test_affine_maps(Float64)
        test_affine_maps(BigFloat)
    end
    @testset "composite maps" begin
        test_composite_maps()
    end
    @testset "product maps" begin
        test_product_map(Float64)
        test_product_map(BigFloat)
    end
    @testset "wrapped maps" begin
        test_wrapped_maps(Float64)
        test_wrapped_maps(BigFloat)
    end
    @testset "scaling maps" begin
        test_scaling_maps(Float64)
        test_scaling_maps(BigFloat)
    end
    @testset "isomorphisms" begin
        test_isomorphisms(Float64)
        test_isomorphisms(BigFloat)
    end
    @testset "Mixed maps" begin
        test_mixed_maps()
    end
end

function test_mixed_maps()
    m1 = composedmap(cos, sin)
    @test domaintype(m1) == Any
    @test m1(0.4) == sin(cos(0.4))

    m2 = composedmap(AffineMap(2.0, 3.0), cos)
    @test domaintype(m2) == Float64
    @test m2(0.4) ≈ cos(2*0.4+3)
    @inferred m2(0.4)
    @test repr(m2) == "cos ∘ (x -> 2.0 * x + 3.0)"

    m3 = composedmap(sin, AffineMap(2.0, 3.0), cos)
    @test m3(0.5) ≈ cos(2*sin(0.5)+3)
    @test repr(m3) == "cos ∘ (x -> 2.0 * x + 3.0) ∘ sin"
    @test domaintype(m3) == Any
    @inferred m3(0.5)

    m4 = productmap(sin, cos)
    @test m4 isa TupleProductMap
    @test domaintype(m4) == Tuple{Any,Any}
    @test m4(0.3,0.5) == (sin(0.3), cos(0.5))
end

test_maps()
