"Supertype of identity maps."
abstract type IdentityMap{T} <: Map{T} end

IdentityMap(n::Int) = DynamicIdentityMap(n)
IdentityMap() = StaticIdentityMap()
IdentityMap(::Val{N}) where {N} = StaticIdentityMap(Val(N))

IdentityMap{T}(n::Int) where {T} = DynamicIdentityMap{T}(n)
IdentityMap{T}(n::Int) where {T<:StaticTypes} = StaticIdentityMap{T}(n)
IdentityMap{T}(::Val{N}) where {N,T} = StaticIdentityMap{T}(Val(N))
IdentityMap{T}() where {T} = StaticIdentityMap{T}()

applymap(map::IdentityMap, x) = x
applymap!(y, map::IdentityMap, x) = y .= x

inverse(m::IdentityMap) = m
inverse(m::IdentityMap, x) = x

islinearmap(::IdentityMap) = true
isrealmap(::IdentityMap{T}) where {T} = isrealtype(T)

isidentitymap(m) = false
isidentitymap(::IdentityMap) = true

mapsize(m::IdentityMap{T}) where {T<:Number} = ()
mapsize(m::IdentityMap{T}) where {T} = (euclideandimension(T),euclideandimension(T))

affinematrix(m::IdentityMap) = identitymatrix(m)
affinevector(m::IdentityMap) = zerovector(m)

jacobian(m::IdentityMap) = ConstantMap(affinematrix(m))
jacobian(m::IdentityMap, x) = affinematrix(m)

jacdet(m::IdentityMap, x) = 1

determinantmap(m::IdentityMap{T}) where {T} = UnityMap{T,prectype(T)}()

mapcompose(m1::IdentityMap) = m1
mapcompose(m1::IdentityMap, maps...) = mapcompose(maps...)
mapcompose2(m1, m2::IdentityMap, maps...) = mapcompose(m1, maps...)

show(io::IO, m::IdentityMap{T}) where {T} = print(io, "x -> x")
map_object_parentheses(m::IdentityMap) = true

"The identity map for variables of type `T`."
struct StaticIdentityMap{T} <: IdentityMap{T}
end

StaticIdentityMap() = StaticIdentityMap{Float64}()
StaticIdentityMap(::Val{N}) where {N} = StaticIdentityMap{SVector{N,Float64}}()

function StaticIdentityMap{T}(n::Int) where {T}
    n == euclideandimension(T) || throw(ArgumentError("Provided dimension of static map is inconsistent with static type"))
    StaticIdentityMap{T}()
end
function StaticIdentityMap{T}(::Val{N}) where {N,T}
    N == euclideandimension(T) || throw(ArgumentError("Provided dimension of static map is inconsistent with static type"))
    StaticIdentityMap{T}()
end

similarmap(m::StaticIdentityMap, ::Type{T}) where {T<:StaticTypes} = StaticIdentityMap{T}()
similarmap(m::StaticIdentityMap, ::Type{T}) where {T} =
    DynamicIdentityMap{T}(euclideandimension(T))

convert(::Type{StaticIdentityMap{T}}, ::StaticIdentityMap) where {T} = StaticIdentityMap{T}()

isequalmap(m1::StaticIdentityMap, m2::StaticIdentityMap) = true
map_hash(m::StaticIdentityMap, h::UInt) = hash("StaticIdentityMap", h)

"Identity map with dynamic size determined by a dimension field."
struct DynamicIdentityMap{T} <: IdentityMap{T}
    dimension   ::  Int
end

const EuclideanIdentityMap{N,T} = StaticIdentityMap{SVector{N,T}}
const VectorIdentityMap{T} = DynamicIdentityMap{Vector{T}}

DynamicIdentityMap(dimension::Int) = VectorIdentityMap(dimension)
VectorIdentityMap(dimension::Int) = VectorIdentityMap{Float64}(dimension)

mapsize(m::DynamicIdentityMap) = (m.dimension, m.dimension)

similarmap(m::DynamicIdentityMap, ::Type{T}) where {T} =
    DynamicIdentityMap{T}(m.dimension)
similarmap(m::DynamicIdentityMap, ::Type{T}) where {T<:StaticTypes} =
    StaticIdentityMap{T}()

isequalmap(m1::DynamicIdentityMap, m2::DynamicIdentityMap) = m1.dimension == m2.dimension
map_hash(m::DynamicIdentityMap, h::UInt) = hashrec("DynamicIdentityMap", m.dimension, h)


"The supertype of constant maps from `T` to `U`."
abstract type ConstantMap{T,U} <: TypedMap{T,U} end

applymap(m::ConstantMap, x) = mapconstant(m)

isconstantmap(m::Map) = false
isconstantmap(m::ConstantMap) = true

isrealmap(m::ConstantMap{T,U}) where {T,U} =
    isrealtype(T) && isrealtype(U) && isreal(mapconstant(m))

mapsize(m::ConstantMap) = _constant_mapsize(m, mapconstant(m))
_constant_mapsize(m::ConstantMap{T,U}, c) where {T<:Number,U<:Number} = ()
_constant_mapsize(m::ConstantMap{T,U}, c) where {T<:Number,U} = (length(c),)
_constant_mapsize(m::ConstantMap{T,U}, c) where {T,U<:Number} = (1,euclideandimension(T))
_constant_mapsize(m::ConstantMap{T,U}, c) where {T,U} = (length(c), euclideandimension(T))

affinematrix(m::ConstantMap) = zeromatrix(m)
affinevector(m::ConstantMap) = mapconstant(m)

jacobian(m::ConstantMap{T}) where {T} = ConstantMap{T}(affinematrix(m))
jacobian(m::ConstantMap, x) = affinematrix(m)

jacdet(::ConstantMap, x) = 0

determinantmap(m::ConstantMap{T}) where {T} = ConstantMap{T}(det(mapconstant(m)))
absmap(m::ConstantMap{T}) where {T} = ConstantMap{T}(abs(mapconstant(m)))

diffvolume(m::ConstantMap{T,U}) where {T,U} = ZeroMap{T,U}()

isequalmap(m1::ConstantMap, m2::ConstantMap) = mapconstant(m1)==mapconstant(m2)
map_hash(m::ConstantMap, h::UInt) = hashrec("ConstantMap", mapconstant(m), h)

similarmap(m::ConstantMap, ::Type{T}) where {T} = ConstantMap{T}(mapconstant(m))
similarmap(m::ConstantMap, ::Type{T}, ::Type{U}) where {T,U} = ConstantMap{T,U}(m.c)

ConstantMap() = ConstantMap{Float64}()
ConstantMap(c) = FixedConstantMap(c)
ConstantMap{T}() where {T} = UnityMap{T}()
ConstantMap{T}(c) where {T} = FixedConstantMap{T}(c)
ConstantMap{T,U}() where {T,U} = UnityMap{T,U}()
ConstantMap{T,U}(c) where {T,U} = FixedConstantMap{T,U}(c)

show(io::IO, m::ConstantMap{T}) where {T} = print(io, "x -> $(mapconstant(m))")
map_object_parentheses(m::ConstantMap) = true


"The zero map `f(x) = 0`."
struct ZeroMap{T,U} <: ConstantMap{T,U}
end
ZeroMap{T}() where {T} = ZeroMap{T,T}()
mapconstant(m::ZeroMap{T,U}) where {T,U} = zero(U)
similarmap(m::ZeroMap{S,U}, ::Type{T}) where {T,S,U} = ZeroMap{T,U}()
similarmap(m::ZeroMap, ::Type{T}, ::Type{U}) where {T,U} = ZeroMap{T,U}()


"The unity map `f(x) = 1`."
struct UnityMap{T,U} <: ConstantMap{T,U}
end
UnityMap{T}() where {T} = UnityMap{T,real(numtype(T))}()
mapconstant(m::UnityMap{T,U}) where {T,U} = one(U)
similarmap(m::UnityMap{S,U}, ::Type{T}) where {T,S,U} = UnityMap{T,U}()
similarmap(m::UnityMap, ::Type{T}, ::Type{U}) where {T,U} = UnityMap{T,U}()


"The constant map `f(x) = c`."
struct FixedConstantMap{T,U} <: ConstantMap{T,U}
    c   ::  U
end
FixedConstantMap{T}(c::U) where {T,U} = FixedConstantMap{T,U}(c)
FixedConstantMap(c::T) where {T} = FixedConstantMap{T}(c)
mapconstant(m::FixedConstantMap) = m.c
