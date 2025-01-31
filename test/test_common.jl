
using FunctionMaps:
    hashrec,
    convert_numtype,
    promote_numtype,
    to_numtype,
    convert_prectype,
    promote_prectype,
    to_prectype,
    convert_eltype,
    euclideandimension,
    isrealtype

function test_hashrec()
    @test hashrec() == 0
    @test hashrec() isa UInt
    A = [1,2]
    @test hashrec(A) == hash(2, hash(1, hash((2,))))
    @test hashrec(1, 2) == hash(1, hash(2))
end

function test_dimension()
    @test euclideandimension(Int) == 1
    @test euclideandimension(Float64) == 1
    @test euclideandimension(ComplexF64) == 1
    @test euclideandimension(SVector{2,Float64}) == 2
    @test euclideandimension(MVector{2,Float64}) == 2
    @test euclideandimension(Tuple{Int,Int}) == 2
    @test euclideandimension(Tuple{Int,Float64}) == 2
    @test_throws ArgumentError euclideandimension(Vector{Float64})
end

function test_prectype()
    @test prectype(:some_symbol) == Any
    @test prectype(1.0) == Float64
    @test prectype(big(1.0)) == BigFloat
    @test prectype(1) == typeof(float(1))
    @test prectype(SVector(1,2)) == typeof(float(1))
    @test prectype(SVector(1,big(2))) == typeof(float(big(2)))
    @test prectype(1.0+2.0im) == Float64
    @test prectype([1.0+2.0im, 3.0]) == Float64
    @test prectype(NTuple{2,Int}) == Float64
    @test prectype(typeof((1.0,))) == Float64
    @test prectype(typeof((1.0, 2.0))) == Float64
    @test prectype(typeof((1.0, 2.0, 3.0))) == Float64
    @test prectype(typeof((1.0, big(2.0), 3.0+im))) == BigFloat
    @test prectype(NTuple{4,Int}) == Float64
    @test @inferred(prectype(1, 2.0)) == Float64
    @test @inferred(prectype(typeof((1, 2.0, 3, 40+im)))) == Float64

    @test convert_prectype(Float64, 2) == 2
    @test convert_prectype(Float64, 2) isa Float64
    @test convert_prectype(BigFloat, 1.0+im) == 1+im
    @test convert_prectype(BigFloat, 1.0+im) isa Complex{BigFloat}
    @test convert_prectype(Float64, SA[1,2]) == SA[1.0,2.0]
    @test convert_prectype(Float64, SA[1,2]) isa SVector{2,Float64}
    @test convert_prectype(Float64, MVector(1,2)) == MVector(1.0,2.0)
    @test convert_prectype(Float64, MVector(1,2)) isa MVector
    @test convert_prectype(Float64, SA[1 2; 3 4]) == SA[1.0 2.0; 3.0 4.0]
    @test convert_prectype(Float64, SA[1 2; 3 4]) isa SMatrix{2,2,Float64}
    @test convert_prectype(Float64, MMatrix{2,2}(1, 2, 3, 4)) == SA[1.0 3.0; 2.0 4.0]
    @test convert_prectype(Float64, MMatrix{2,2}(1, 2, 3, 4)) isa MMatrix{2,2,Float64}
    @test convert_prectype(BigFloat, SA[1,2+im]) isa SVector{2,Complex{BigFloat}}
    @test_throws ArgumentError convert_prectype(BigFloat, "a")

    @test promote_prectype(2) == 2
    @test promote_prectype(2, 3.0) isa Tuple{Float64,Float64}
    @test promote_prectype(2, 3.0+im, big(4)) isa Tuple{BigFloat,Complex{BigFloat},BigFloat}

    @test to_prectype(Float64, Int) === Float64
    @test to_prectype(Float32, Float64) === Float32
    @test to_prectype(Int, Int) === Int
    
    @test to_prectype(Float64, Complex{Int}) === Complex{Float64}
    @test to_prectype(Float32, Complex{Float64}) === Complex{Float32}
    @test to_prectype(Float64, Vector{Int}) === Vector{Float64}
    @test to_prectype(Float32, Vector{Float64}) === Vector{Float32}
    @test to_prectype(Float64, SVector{3,Int}) === SVector{3,Float64}
    @test to_prectype(Float32, MVector{3,Float64}) === MVector{3,Float32}
    @test to_prectype(Float64, SMatrix{2,2,Int}) === SMatrix{2,2,Float64}
    @test to_prectype(Float32, MMatrix{2,2,Float64}) === MMatrix{2,2,Float32}
    @test_throws ArgumentError to_prectype(Float64, String)
    @test_throws ArgumentError to_prectype(Float64, Dict{Int,Float64})    
end


function test_numtype()
    @test numtype(:some_symbol) == Any
    @test numtype(1.0) == Float64
    @test numtype(big(1.0)) == BigFloat
    @test numtype(1) == Int
    @test numtype([1,2,3], [4,5,6]) == Int
    @test numtype(SVector(1,2)) == Int
    @test numtype(SVector(1,big(2))) == BigInt
    @test numtype(Array{Float64,2}) == Float64
    @test numtype(1.0+2.0im) == Complex{Float64}
    @test numtype([1.0+2.0im, 3.0]) == Complex{Float64}
    @test numtype(NTuple{2,Complex{Int}}) == Complex{Int}
    @test numtype(typeof((1.0,))) == Float64
    @test numtype(typeof((1.0, 2.0))) == Float64
    @test numtype((1.0, 2.0, 3.0)) == Float64
    @test numtype((1.0, 2.0, 3.0, 4.0)) == Float64
    @test numtype(1.0, big(2.0), 3.0+im) == Complex{BigFloat}
    @test numtype(typeof((1.0, big(2.0), 3.0+im))) == Complex{BigFloat}
    @test numtype(Tuple{Int,Int,Int,Int,Float64}) == Float64
    @test @inferred(numtype(1, 2.0)) == Float64
    @test @inferred(numtype(typeof((1, 2.0, 3, 40+im)))) == Complex{Float64}

    @test convert_numtype(Float64, 2) == 2
    @test convert_numtype(Float64, 2) isa Float64
    @test convert_numtype(Float64, SA[1,2]) == SA[1.0,2.0]
    @test convert_numtype(Float64, SA[1,2]) isa SVector{2,Float64}
    @test convert_numtype(Float64, MVector(1,2)) == MVector(1.0,2.0)
    @test convert_numtype(Float64, MVector(1,2)) isa MVector
    @test convert_numtype(Float64, SA[1 2; 3 4]) == SA[1.0 2.0; 3.0 4.0]
    @test convert_numtype(Float64, SA[1 2; 3 4]) isa SMatrix{2,2,Float64}
    @test convert_numtype(Float64, MMatrix{2,2}(1, 2, 3, 4)) == SA[1.0 3.0; 2.0 4.0]
    @test convert_numtype(Float64, MMatrix{2,2}(1, 2, 3, 4)) isa MMatrix{2,2,Float64}
    @test_throws ArgumentError convert_numtype(BigFloat, "a")

    @test promote_numtype(2) == 2
    @test promote_numtype(2, 3.0) isa Tuple{Float64,Float64}
    @test promote_numtype(2, 3.0+im, big(4)) isa Tuple{Complex{BigFloat},Complex{BigFloat},Complex{BigFloat}}

    @test to_numtype(Float64, Int) === Float64
    @test to_numtype(Float32, Float64) === Float32
    @test to_numtype(Int, Int) === Int
    
    @test to_numtype(Float64, Vector{Int}) === Vector{Float64}
    @test to_numtype(Float32, Vector{Float64}) === Vector{Float32}
    @test to_numtype(Float64, SVector{3,Int}) === SVector{3,Float64}
    @test to_numtype(Float32, MVector{3,Float64}) === MVector{3,Float32}
    @test to_numtype(Float64, SMatrix{2,2,Int}) === SMatrix{2,2,Float64}
    @test to_numtype(Float32, MMatrix{2,2,Float64}) === MMatrix{2,2,Float32}
    @test_throws ArgumentError to_numtype(Float64, String)
    @test_throws ArgumentError to_numtype(Float64, Dict{Int,Float64})
end

function test_eltype()
    @test convert_eltype(Float64, Diagonal([1,2])) == Diagonal([1,2])
    @test eltype(convert_eltype(Float64, Diagonal([1,2]))) == Float64
    @test convert_eltype(Float64, 1:4) == 1:4
    @test eltype(convert_eltype(Float64, 1:4)) == Float64
    @test convert_eltype(Float64, Set([1,4])) == Set([1,4])
    @test eltype(convert_eltype(Float64, Set([1,4]))) == Float64
    @test convert_eltype(Float64, 5) == 5
    @test eltype(convert_eltype(Float64, 5)) == Float64

    @test FunctionMaps.promotable_eltypes(Int,Float64)
    @test FunctionMaps.promotable_eltypes(Vector{Int},Vector{Float64})
end

function test_realtype()
    @test isrealtype(Any) == false
    @test isrealtype(Int) == true
    @test isrealtype(ComplexF64) == false
    @test isrealtype(Matrix{Float64}) == true
end

@testset "common functionality" begin
    test_hashrec()
    test_dimension()
    test_prectype()
    test_numtype()
    test_eltype()
    test_realtype()
end
