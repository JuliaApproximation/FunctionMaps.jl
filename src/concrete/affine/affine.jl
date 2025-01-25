"""
    AffineMap{T} <: AbstractAffineMap{T}

    The supertype of all affine maps that store `A` and `b`.
Concrete subtypes differ in how `A` and `b` are represented.
"""
abstract type AffineMap{T} <: AbstractAffineMap{T} end

"""
    AffineMap(A, b)

Return an affine map with an appropriate concrete type depending on the arguments
`A` and `b`.

# Examples
```julia
julia> AffineMap(2, 3)
x -> 2 * x + 3
```
"""
AffineMap(A::Number, b::Number) = ScalarAffineMap(A, b)
AffineMap(A::StaticMatrix, b::StaticVector) = StaticAffineMap(A, b)
AffineMap(A::Matrix, b::Vector) = VectorAffineMap(A, b)
AffineMap(A::UniformScaling{Bool}, b::Number) = ScalarAffineMap(one(b), b)
AffineMap(A, b) = GenericAffineMap(A, b)

AffineMap{T}(A::Number, b::Number) where {T<:Number} = ScalarAffineMap{T}(A, b)
AffineMap{T}(A::AbstractMatrix, b::AbstractVector) where {N,S,T<:SVector{N,S}} = StaticAffineMap{S,N}(A, b)
AffineMap{T}(A::Matrix, b::Vector) where {S,T<:Vector{S}} = VectorAffineMap{S}(A, b)
AffineMap{T}(A::UniformScaling{Bool}, b::Number) where {T} = ScalarAffineMap{T}(one(T), b)
AffineMap{T}(A, b) where {T} = GenericAffineMap{T}(A, b)

similarmap(m::AffineMap, ::Type{T}) where {T} = AffineMap{T}(m.A, m.b)

convert(::Type{AffineMap}, m) = (@assert isaffinemap(m); AffineMap(affinematrix(m), affinevector(m)))
convert(::Type{AffineMap{T}}, m) where {T} = (@assert isaffinemap(m); AffineMap{T}(affinematrix(m), affinevector(m)))
# avoid ambiguity errors with convert(::Type{T}, x::T) in Base:
convert(::Type{AffineMap}, m::AffineMap) = m
convert(::Type{AffineMap{T}}, m::AffineMap{T}) where T = m

# If y = A*x+b, then x = inv(A)*(y-b) = inv(A)*y - inv(A)*b
inverse(m::AffineMap) = (@assert issquaremap(m); AffineMap(inv(m.A), -inv(m.A)*m.b))
inverse(m::AffineMap, x) = (@assert issquaremap(m); m.A \ (x-m.b))

function leftinverse(m::AffineMap)
    @assert isoverdetermined(m)
    pA = matrix_pinv(m.A)
    AffineMap(pA, -pA*m.b)
end
function rightinverse(m::AffineMap)
    @assert isunderdetermined(m)
    pA = matrix_pinv(m.A)
    AffineMap(pA, -pA*m.b)
end
function leftinverse(m::AffineMap, x)
    @assert isoverdetermined(m)
    m.A \ (x-m.b)
end
function rightinverse(m::AffineMap, x)
    @assert isunderdetermined(m)
    m.A \ (x-m.b)
end


"An affine map for any combination of types of `A` and `b`."
struct GenericAffineMap{T,AA,B} <: AffineMap{T}
    A   ::  AA
    b   ::  B
end

GenericAffineMap(A, b) = GenericAffineMap{typeof(b)}(A, b)
GenericAffineMap(A::AbstractVector{S}, b::AbstractVector{T}) where {S,T} =
    GenericAffineMap{promote_type(S,T)}(A, b)
GenericAffineMap(A::AbstractArray{S}, b::AbstractVector{T}) where {S,T} =
    GenericAffineMap{Vector{promote_type(S,T)}}(A, b)
GenericAffineMap(A::StaticMatrix{M,N,S}, b::StaticVector{M,T}) where {M,N,S,T} =
    GenericAffineMap{SVector{N,promote_type(S,T)}}(A, b)
GenericAffineMap(A::StaticMatrix{M,N,S}, b::AbstractVector{T}) where {M,N,S,T} =
    GenericAffineMap{SVector{N,promote_type(S,T)}}(A, b)
GenericAffineMap(A::S, b::AbstractVector{T}) where {S<:Number,T} =
    GenericAffineMap{Vector{promote_type(S,T)}}(A, b)
GenericAffineMap(A::S, b::StaticVector{N,T}) where {S<:Number,N,T} =
    GenericAffineMap{SVector{N,promote_type(S,T)}}(A, b)
GenericAffineMap(A::UniformScaling{Bool}, b) =
    GenericAffineMap(UniformScaling{eltype(b)}(1), b)


# Fallback routine for generic A and b, special cases follow
GenericAffineMap{T}(A, b) where {T} = GenericAffineMap{T,typeof(A),typeof(b)}(A, b)

GenericAffineMap{T}(A::AbstractVector{S}, b::AbstractVector{U}) where {T<:Number,S,U} =
    GenericAffineMap{T}(convert(AbstractVector{T}, A), convert(AbstractVector{T}, b))
GenericAffineMap{T}(A::AbstractVector{T}, b::AbstractVector{T}) where {T<:Number} =
    GenericAffineMap{T,typeof(A),typeof(b)}(A, b)
GenericAffineMap{T}(A::Number, b) where {T} = GenericAffineMap{T,eltype(T),typeof(b)}(A, b)
GenericAffineMap{T}(A::Number, b::AbstractVector) where {N,S,T <: StaticVector{N,S}} =
    GenericAffineMap{T,S,SVector{N,S}}(A, b)
# Promote element types of abstract arrays
GenericAffineMap{T}(A::AbstractArray, b::AbstractVector) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T}(convert(AbstractArray{eltype(T)},A), convert(AbstractVector{eltype(T)}, b))
GenericAffineMap{T}(A::AbstractArray{S}, b::AbstractVector{S}) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T,typeof(A),typeof(b)}(A, b)
GenericAffineMap{T}(A::UniformScaling{Bool}, b::AbstractVector) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T}(A*one(S), convert(AbstractVector{S}, b))
GenericAffineMap{T}(A::UniformScaling{S}, b::AbstractVector{S}) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T,typeof(A),typeof(b)}(A, b)


similarmap(m::GenericAffineMap, ::Type{T}) where {T} = AffineMap{T}(m.A, m.b)

convert(::Type{GenericAffineMap{T}}, m::GenericAffineMap) where {T} =
    GenericAffineMap{T}(m.A, m.b)



"An affine map with scalar representation."
struct ScalarAffineMap{T} <: AffineMap{T}
    A   ::  T
    b   ::  T
end

ScalarAffineMap(A, b) = ScalarAffineMap(promote(A, b)...)

isrealmap(m::ScalarAffineMap{T}) where {T} = isrealtype(T)

show(io::IO, m::ScalarAffineMap) = show_scalar_affine_map(io, m.A, m.b)
show_scalar_affine_map(io, A::Real, b::Real) = print(io, "x -> $(A) * x", b < 0 ? " - " : " + ", abs(b))
show_scalar_affine_map(io, A::Complex, b::Complex) = print(io, "x -> ($(A)) * x + ", b)
show_scalar_affine_map(io, A, b) = print(io, "x -> ($(A)) * x + $(b)")


convert(::Type{ScalarAffineMap{T}}, m::ScalarAffineMap) where {T} =
    ScalarAffineMap{T}(m.A, m.b)

"An affine map with array and vector representation."
struct VectorAffineMap{T} <: AffineMap{Vector{T}}
    A   ::  Matrix{T}
    b   ::  Vector{T}
end

VectorAffineMap(A::AbstractArray{T}, b::AbstractVector{T}) where {T} =
    VectorAffineMap{T}(A, b)
function VectorAffineMap(A::AbstractArray{S}, b::AbstractVector{T}) where {S,T}
    U = promote_type(S,T)
    VectorAffineMap(convert(AbstractArray{U}, A), convert(AbstractVector{U}, b))
end

convert(::Type{VectorAffineMap{T}}, m::VectorAffineMap) where {T} =
    VectorAffineMap{T}(m.A, m.b)



"An affine map with representation using static arrays."
struct StaticAffineMap{T,N,M,L} <: AffineMap{SVector{N,T}}
    A   ::  SMatrix{M,N,T,L}
    b   ::  SVector{M,T}
end

# Constructors:
# - first, we deduce T
StaticAffineMap(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T} =
    StaticAffineMap{T}(A, b)
function StaticAffineMap(A::AbstractMatrix{S}, b::AbstractVector{T}) where {S,T}
    U = promote_type(S,T)
    StaticAffineMap(convert(AbstractMatrix{U}, A), convert(AbstractVector{U}, b))
end

StaticAffineMap{T}(A::AbstractMatrix, b::AbstractVector) where {T} =
    StaticAffineMap{T}(convert(AbstractMatrix{T}, A), convert(AbstractVector{T}, b))

# - then, we determine N and/or M, from the arguments
function StaticAffineMap{T}(A::AbstractMatrix{T}, b::StaticVector{M,T}) where {T,M}
    @assert size(A) == (M,M)
    StaticAffineMap{T,M,M}(A, b)
end
StaticAffineMap{T}(A::StaticMatrix{M,N,T}, b::AbstractVector) where {T,N,M} =
    StaticAffineMap{T,N,M}(A, b)
StaticAffineMap{T}(A::StaticMatrix{M,N,T}, b::StaticVector{M,T}) where {T,N,M} =
    StaticAffineMap{T,N,M}(A, b)
# line below catches ambiguity error
StaticAffineMap{T}(A::StaticMatrix{M1,N,T}, b::StaticVector{M2,T}) where {T,N,M1,M2} =
    throw(ArgumentError("Non-matching dimensions"))
StaticAffineMap{T,N}(A::AbstractMatrix, b::AbstractVector) where {T,N} =
    StaticAffineMap{T,N,N}(A, b)
StaticAffineMap{T,N}(A::StaticMatrix{M,N}, b::AbstractVector) where {T,N,M} =
    StaticAffineMap{T,N,M}(A, b)

# - finally invoke the constructor (and implicitly convert the data if necessary)
StaticAffineMap{T,N,M}(A::AbstractMatrix, b::AbstractVector) where {T,N,M} =
    StaticAffineMap{T,N,M,M*N}(A, b)

convert(::Type{Map{SVector{N,T}}}, m::VectorAffineMap) where {N,T} =
    StaticAffineMap{T,N}(m.A, m.b)
