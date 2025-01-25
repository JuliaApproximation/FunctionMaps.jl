
"A `Translation` represents the map `y = x + b`."
abstract type Translation{T} <: AbstractAffineMap{T} end

# unsafe_matrix(m::Translation) = I
affinematrix(m::Translation) = to_matrix(domaintype(m), LinearAlgebra.I, unsafe_vector(m))
affinevector(m::Translation) = unsafe_vector(m)

mapsize(m::Translation) = _translation_mapsize(m, domaintype(m), unsafe_vector(m))
_translation_mapsize(m, ::Type{T}, b::Number) where {T} = ()
_translation_mapsize(m, ::Type{T}, b) where {T} = (length(b),length(b))

map_stencil(m::Translation, x) = _translation_map_stencil(m, x, unsafe_vector(m))
_translation_map_stencil(m, x, b) = [x, " + ", b]
_translation_map_stencil(m, x, b::Real) =
    b >= 0 ? [x, " + ", b] : [x, " - ", abs(b)]
map_stencil_broadcast(m::Translation, x) = _translation_map_stencil_broadcast(m, x, unsafe_vector(m))
_translation_map_stencil_broadcast(m, x, b) = [x, " .+ ", b]

"Translation by a scalar value."
struct ScalarTranslation{T} <: Translation{T}
    b   ::  T
end

isrealmap(m::ScalarTranslation{T}) where {T} = isrealtype(T)

show(io::IO, m::ScalarTranslation) = show_scalar_translation(io, m.b)
show_scalar_translation(io, b::Real) = print(io, "x -> x", b < 0 ? " - " : " + ", abs(b))
show_scalar_translation(io, b) = print(io, "x -> x + ", b)


"Translation by a static vector."
struct StaticTranslation{T,N} <: Translation{SVector{N,T}}
    b   ::  SVector{N,T}
end

"Translation by a vector."
struct VectorTranslation{T} <: Translation{Vector{T}}
    b   ::  Vector{T}
end

"Translation by a generic vectorlike object."
struct GenericTranslation{T,B} <: Translation{T}
    b   ::  B
end

Translation(b::Number) = ScalarTranslation(b)
Translation(b::StaticVector) = StaticTranslation(b)
Translation(b::Vector) = VectorTranslation(b)
Translation(b) = GenericTranslation(b)

Translation{T}(b::Number) where {T<:Number} = ScalarTranslation{T}(b)
Translation{T}(b::AbstractVector) where {N,S,T<:StaticVector{N,S}} = StaticTranslation{S,N}(b)
Translation{T}(b::Vector) where {S,T<:Vector{S}} = VectorTranslation{S}(b)
Translation{T}(b) where {T} = GenericTranslation{T}(b)

jacdet(m::Translation{T}) where {T} = UnityMap{T,eltype(T)}()
jacdet(m::Translation{T}, x) where {T} = one(eltype(T))

isrealmap(m::Translation) = isreal(unsafe_vector(m))

isequalmap(m1::Translation, m2::Translation) = unsafe_vector(m1)==unsafe_vector(m2)

similarmap(m::Translation, ::Type{T}) where {T} = Translation{T}(m.b)

applymap(m::Translation, x) = _translation_applymap(m, x, unsafe_vector(m))
_translation_applymap(m, x, b) = x + b
applymap!(y, m::Translation, x) = _translation_applymap!(y, m, x, unsafe_vector(m))
_translation_applymap!(y, m, x, b) = y .= x .+ m.b

inverse(m::Translation{T}) where {T} = Translation{T}(-m.b)
inverse(m::Translation, x) = x - m.b

ScalarTranslation(b::Number) = ScalarTranslation{typeof(b)}(b)

StaticTranslation(b::AbstractVector{T}) where {T} = StaticTranslation{T}(b)

StaticTranslation{T}(b::StaticVector{N}) where {N,T} =
    StaticTranslation{T,N}(b)

VectorTranslation(b::AbstractVector{T}) where {T} = VectorTranslation{T}(b)

GenericTranslation(b) = GenericTranslation{typeof(b)}(b)
GenericTranslation(b::AbstractVector{T}) where {T} =
    GenericTranslation{Vector{T}}(b)
GenericTranslation(b::StaticVector{N,T}) where {N,T} =
    GenericTranslation{SVector{N,T}}(b)

GenericTranslation{T}(b) where {T} = GenericTranslation{T,typeof(b)}(b)
GenericTranslation{T}(b::Number) where {T<:Number} =
    GenericTranslation{T,T}(b)
