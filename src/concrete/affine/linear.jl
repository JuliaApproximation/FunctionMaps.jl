"""
    LinearMap{T} <: AbstractAffineMap{T}

The supertype of all linear maps `y = A*x`.

Concrete subtypes may differ in how `A` is represented.
"""
abstract type LinearMap{T} <: AbstractAffineMap{T} end

mapsize(m::LinearMap) = _linearmap_size(m, domaintype(m), unsafe_matrix(m))
_linearmap_size(m, T, A) = size(A)
_linearmap_size(m, ::Type{T}, A::Number) where {T<:Number} = ()
_linearmap_size(m, ::Type{T}, A::AbstractVector) where {T<:Number} = (length(A),)
_linearmap_size(m, ::Type{T}, A::Number) where {N,T<:StaticVector{N}} = (N,N)

affinematrix(m::LinearMap) = to_matrix(domaintype(m), unsafe_matrix(m))
affinevector(m::LinearMap) = to_vector(domaintype(m), unsafe_matrix(m))

LinearMap(A::Number) = ScalarLinearMap(A)
LinearMap(A::SMatrix) = StaticLinearMap(A)
LinearMap(A::MMatrix) = StaticLinearMap(A)
LinearMap(A::Matrix) = VectorLinearMap(A)
LinearMap(A) = GenericLinearMap(A)

LinearMap{T}(A::Number) where {T<:Number} = _LinearMap(A, T, promote_type(T,typeof(A)))
_LinearMap(A::Number, ::Type{T}, ::Type{T}) where {T<:Number} = ScalarLinearMap{T}(A)
_LinearMap(A::Number, ::Type{T}, ::Type{S}) where {T<:Number,S} = GenericLinearMap{T}(A)
LinearMap{T}(A::SMatrix{M,N}) where {M,N,S,T <: SVector{N,S}} = StaticLinearMap{S}(A)
LinearMap{T}(A::MMatrix{M,N}) where {M,N,S,T <: SVector{N,S}} = StaticLinearMap{S}(A)
LinearMap{T}(A::Matrix) where {S,T <: Vector{S}} = VectorLinearMap{S}(A)
LinearMap{T}(A) where {T} = GenericLinearMap{T}(A)

# convenience functions
LinearMap(a::Number...) = LinearMap(promote(a...))
LinearMap(a::NTuple{N,T}) where {N,T <: Number} = LinearMap{SVector{N,T}}(Diagonal(SVector{N,T}(a)))

applymap(m::LinearMap, x) = _linear_applymap(m, x, unsafe_matrix(m))
_linear_applymap(m, x, A) = A*x

applymap!(y, m::LinearMap, x) = _linear_applymap!(y, m, x, unsafe_matrix(m))
_linear_applymap!(y, m, x, A) = mul!(y, A, x)

islinearmap(m::LinearMap) = true

isrealmap(m::LinearMap) = _linear_isrealmap(m, unsafe_matrix(m))
_linear_isrealmap(m, A) = isreal(A)

isequalmap(m1::LinearMap, m2::LinearMap) = affinematrix(m1) == affinematrix(m2)

# inverse should be called only on square maps, otherwise use
# leftinverse or rightinverse in order to use pinv instead of inv
inverse(m::LinearMap) = (@assert issquaremap(m); LinearMap(inv(m.A)))
inverse(m::LinearMap, x) = (@assert issquaremap(m); m.A \ x)

function leftinverse(m::LinearMap)
    @assert isoverdetermined(m)
    LinearMap(matrix_pinv(m.A))
end
function rightinverse(m::LinearMap)
    @assert isunderdetermined(m)
    LinearMap(matrix_pinv(m.A))
end
function leftinverse(m::LinearMap, x)
    @assert isoverdetermined(m)
    m.A \ x
end
function rightinverse(m::LinearMap, x)
    @assert isunderdetermined(m)
    m.A \ x
end

similarmap(m::LinearMap, ::Type{T}) where {T} = LinearMap{T}(m.A)

convert(::Type{Map}, a::Number) = LinearMap(a)
convert(::Type{Map{T}}, a::Number) where {T} = LinearMap{T}(a)

convert(::Type{LinearMap}, ::StaticIdentityMap{T}) where {T} = LinearMap{T}(1)
convert(::Type{LinearMap{T}}, ::StaticIdentityMap) where {T} = LinearMap{T}(1)
convert(::Type{AbstractAffineMap}, m::StaticIdentityMap) = convert(LinearMap, m)
convert(::Type{AbstractAffineMap{T}}, m::StaticIdentityMap) where {T} = convert(LinearMap, m)


map_stencil(m::LinearMap, x) = [unsafe_matrix(m), " * ", x]
map_stencil_broadcast(m::LinearMap, x) = _linear_map_stencil_broadcast(m, x, unsafe_matrix(m))
_linear_map_stencil_broadcast(m, x, A) = [A, " .* ", x]
_linear_map_stencil_broadcast(m, x, A::Number) = [A, " * ", x]


"A `GenericLinearMap` is a linear map `y = A*x` for any type of `A`."
struct GenericLinearMap{T,AA} <: LinearMap{T}
    A   ::  AA
end

"""
What is the suggested domaintype for a generic linear map `A*x` with
the given argument 'A'?
"""
glm_domaintype(A) = Any
glm_domaintype(A::Number) = typeof(A)
glm_domaintype(A::AbstractMatrix{T}) where T = Vector{T}
glm_domaintype(A::StaticMatrix{M,N,T}) where {M,N,T} = SVector{N,T}
glm_domaintype(A::AbstractVector{T}) where {T} = T
glm_domaintype(A::Diagonal{T,<:StaticVector{N,T}}) where {N,T} = SVector{N,T}

GenericLinearMap(A) = GenericLinearMap{glm_domaintype(A)}(A)

# Allow any A
GenericLinearMap{T}(A) where {T} = GenericLinearMap{T,typeof(A)}(A)

# Promote some eltypes if applicable
GenericLinearMap{T}(A::Number) where {T <: Number} =
    _GenericLinearMap(A, T, promote_type(T,typeof(A)))
_GenericLinearMap(A::Number, ::Type{T}, ::Type{T}) where {T <: Number} =
    GenericLinearMap{T,T}(A)
_GenericLinearMap(A::Number, ::Type{T}, ::Type{S}) where {S,T<:Number} =
    GenericLinearMap{T,typeof(A)}(A)

GenericLinearMap{T}(A::Number) where {T <: AbstractVector} =
    _GenericLinearMap(A, T, promote_type(eltype(T),typeof(A)))
_GenericLinearMap(A::Number, ::Type{T}, ::Type{S}) where {S,T<:AbstractVector{S}} =
    GenericLinearMap{T,S}(A)
_GenericLinearMap(A::Number, ::Type{T}, ::Type{S}) where {S,T<:AbstractVector} =
    GenericLinearMap{T,typeof(A)}(A)

GenericLinearMap{T}(A::AbstractMatrix{S}) where {S,T <: AbstractVector{S}} =
    GenericLinearMap{T,typeof(A)}(A)
GenericLinearMap{T}(A::AbstractMatrix{U}) where {S,T <: AbstractVector{S},U} =
    GenericLinearMap{T}(convert(AbstractMatrix{S}, A))
GenericLinearMap{T}(A::AbstractVector{T}) where {T} =
    GenericLinearMap{T,typeof(A)}(A)
GenericLinearMap{T}(A::AbstractVector{S}) where {S,T} =
    GenericLinearMap{T}(convert(AbstractVector{T}, A))
GenericLinearMap{T}(A::AbstractMatrix) where {T<:Number} =
    throw(ArgumentError("Linear map with matrix A can not have scalar domaintype."))

# Preserve the action on vectors with a number type
inverse(m::GenericLinearMap{T,AA}) where {T<:AbstractVector,AA<:Number} =
    LinearMap{T}(inv(m.A))
leftinverse(m::GenericLinearMap{T,AA}) where {T<:AbstractVector,AA<:Number} =
    LinearMap{T}(inv(m.A))
rightinverse(m::GenericLinearMap{T,AA}) where {T<:AbstractVector,AA<:Number} =
    LinearMap{T}(inv(m.A))

"A `ScalarLinearMap` is a linear map `y = A*x` for scalars."
struct ScalarLinearMap{T} <: LinearMap{T}
    A   ::  T
end

ScalarLinearMap(A::Number) = ScalarLinearMap{typeof(A)}(A)

isrealmap(m::ScalarLinearMap{T}) where {T} = isrealtype(T)

show(io::IO, m::ScalarLinearMap) = show_scalar_linear_map(io, m.A)
show_scalar_linear_map(io, A::Real) = print(io, "x -> $(A) * x")
show_scalar_linear_map(io, A::Complex) = print(io, "x -> ($(A)) * x")


"A `VectorLinearMap` is a linear map `y = A*x` using vectors and matrices."
struct VectorLinearMap{T} <: LinearMap{Vector{T}}
    A   ::  Matrix{T}
end

VectorLinearMap(A::AbstractMatrix{T}) where {T} =
    VectorLinearMap{T}(A)


"A `StaticLinearMap` is a linear map `y = A*x` using static arrays."
struct StaticLinearMap{T,N,M,L} <: LinearMap{SVector{N,T}}
    A   ::  SMatrix{M,N,T,L}
end

StaticLinearMap(A::AbstractMatrix{T}) where {T} =
    StaticLinearMap{T}(A)

StaticLinearMap{T}(A::StaticMatrix{M,N,S}) where {M,N,T,S} =
    StaticLinearMap{T}(convert(AbstractMatrix{T}, A))
StaticLinearMap{T}(A::StaticMatrix{M,N,T}) where {M,N,T} =
    StaticLinearMap{T,M,N}(A)
StaticLinearMap{T,N,M}(A::AbstractMatrix) where {T,N,M} =
    StaticLinearMap{T,N,M,M*N}(A)

convert(::Type{Map{SVector{N,T}}}, m::VectorLinearMap) where {N,T} = StaticLinearMap{T,N,N}(m.A)



