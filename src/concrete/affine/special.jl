
##############################
# Numbers and arrays as a map
##############################

# Arrays represent the linear map A*x
MapStyle(::Type{<:AbstractArray}) = IsMap()

Map(A::AbstractArray) = LinearMap(A)
Map{T}(A::AbstractArray) where T = LinearMap{T}(A)

domaintype(A::AbstractArray) = domaintype(Map(A))

applymap(A::AbstractArray, x) = A*x
mapsize(A::AbstractArray) = size(A)

islinearmap(A::AbstractArray) = true
isaffinemap(A::AbstractArray) = true
affinematrix(A::AbstractArray) = A
affinevector(A::AbstractArray) = zerovector(A)

inverse(A::AbstractMatrix) = inv(A)
inverse(A::AbstractMatrix, x) = A \ x

jacobian(A::AbstractMatrix) = ConstantMap{glm_domaintype(A)}(A)
jacobian(A::AbstractMatrix, x) = A

canonicalmap(A::AbstractArray) = Map(A)
canonicalmap(::Equal, A::AbstractArray) = Map(A)

# Numbers represent the linear map a*x

MapStyle(::Type{<:Number}) = IsMap()

applymap(a::Number, x) = a*x

Map(a::Number) = LinearMap(a)
Map{T}(a::Number) where T = LinearMap{T}(a)
domaintype(a::Number) = typeof(a)

mapsize(a::Number) = ()

islinearmap(a::Number) = true
isaffinemap(a::Number) = true
affinematrix(a::Number) = a
affinevector(a::Number) = zero(a)

inverse(a::Number) = inv(a)
inverse(a::Number, x) = a \ x

jacobian(a::Number) = ConstantMap(a)
jacobian(a::Number, x) = a

canonicalmap(a::Number) = Map(a)
canonicalmap(::Equal, a::Number) = Map(a)


###########################
# Uniform scaling as a map
###########################

MapStyle(::Type{<:UniformScaling}) = IsMap()

convert(::Type{Map}, A::UniformScaling) = GenericLinearMap{Vector{Any}}(A)
convert(::Type{Map{T}}, A::UniformScaling) where {T} = GenericLinearMap{T}(A)

# I(4) returns a diagonal matrix: the action of I is multiplication
applymap(I::UniformScaling, x) = I*x
domaintype(::UniformScaling{T}) where T = T
islinearmap(::UniformScaling) = true
affinematrix(I::UniformScaling) = I
affinevector(I::UniformScaling) = zerovector(I)
