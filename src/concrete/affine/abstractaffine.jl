
"""
    AbstractAffineMap{T} <: Map{T}

An affine map has the general form `y = A*x + b`.

We use `affinematrix(m)` and `affinevector(m)` to denote `A` and `b`
respectively. Concrete subtypes include linear maps of the form `y = A*x`
and translations of the form `y = x + b`.

See also: [`affinematrix`](@ref), [`affinevector`](@ref).
"""
abstract type AbstractAffineMap{T} <: Map{T} end

unsafe_matrix(m::AbstractAffineMap) = m.A
unsafe_vector(m::AbstractAffineMap) = m.b

"Return the matrix `A` in the affine map `Ax+b`."
affinematrix(m::AbstractAffineMap) = to_matrix(domaintype(m), unsafe_matrix(m), unsafe_vector(m))

"Return the vector `b` in the affine map `Ax+b`."
affinevector(m::AbstractAffineMap) = to_vector(domaintype(m), unsafe_matrix(m), unsafe_vector(m))

applymap(m::AbstractAffineMap, x) = _affine_applymap(m, x, unsafe_matrix(m), unsafe_vector(m))
_affine_applymap(m, x, A, b) = A*x + b

applymap!(y, m::AbstractAffineMap, x) = _affine_applymap!(y, m, x, unsafe_matrix(m), unsafe_vector(m))
function _affine_applymap!(y, m, x, A, b)
    mul!(y, A, x)
    y .+= b
    y
end

isrealmap(m::AbstractAffineMap) = _affine_isrealmap(m, unsafe_matrix(m), unsafe_vector(m))
_affine_isrealmap(m, A, b) = isrealmap(A) && isreal(b)

jacobian(m::AbstractAffineMap{T}) where {T} = ConstantMap{T}(affinematrix(m))
jacobian(m::AbstractAffineMap, x) = affinematrix(m)

jacdet(m::AbstractAffineMap, x) = _affine_jacdet(m, x, unsafe_matrix(m))
_affine_jacdet(m, x, A) = det(A)
_affine_jacdet(m, x::Number, A::UniformScaling) = A.λ
_affine_jacdet(m, x::AbstractVector, A::Number) = A^length(x)
_affine_jacdet(m, x::AbstractVector, A::UniformScaling) = A.λ^length(x)

function diffvolume(m::AbstractAffineMap{T}) where T
    J = jacobian(m)
    c = sqrt(det(affinevector(J)'*affinevector(J)))
    ConstantMap{T}(c)
end

"""
    islinearmap(m)

Is `m` a linear map?
"""
islinearmap(m) = false
islinearmap(m::AbstractAffineMap) = _affine_islinearmap(m, unsafe_vector(m))
_affine_islinearmap(m, b) = all(b .== 0)

"""
    isaffinemap(m)

Is `m` an affine map?

If `m` is affine, then it has the form `m(x) = A*x+b`.

See also: [`affinematrix`](@ref), [`affinevector`](@ref).
"""
isaffinemap(m) = false
isaffinemap(m::Map) = islinearmap(m) || isconstantmap(m)
isaffinemap(m::AbstractAffineMap) = true

isequalmap(m1::AbstractAffineMap, m2::AbstractAffineMap) =
    affinematrix(m1) == affinematrix(m2) && affinevector(m1) == affinevector(m2)

isequalmap(m1::AbstractAffineMap, m2::IdentityMap) =
    islinearmap(m1) && affinematrix(m1) == affinematrix(m2)
isequalmap(m1::IdentityMap, m2::AbstractAffineMap) = isequalmap(m2, m1)

map_hash(m::AbstractAffineMap, h::UInt) = hashrec("AbstractAffineMap", affinematrix(m), affinevector(m), h)

mapsize(m::AbstractAffineMap) = _affine_mapsize(m, domaintype(m), unsafe_matrix(m), unsafe_vector(m))
_affine_mapsize(m, T, A::AbstractArray, b) = size(A)
_affine_mapsize(m, T, A::AbstractVector, b::AbstractVector) = (length(A),)
_affine_mapsize(m, T, A::Number, b::Number) = ()
_affine_mapsize(m, T, A::Number, b::AbstractVector) = (length(b),length(b))
_affine_mapsize(m, T, A::UniformScaling, b::Number) = ()
_affine_mapsize(m, T, A::UniformScaling, b) = (length(b),length(b))


Display.displaystencil(m::AbstractAffineMap) = vcat(["x -> "], map_stencil(m, 'x'))
show(io::IO, mime::MIME"text/plain", m::AbstractAffineMap) = composite_show(io, mime, m)

map_stencil(m::AbstractAffineMap, x) = _affine_map_stencil(m, x, unsafe_matrix(m), unsafe_vector(m))
_affine_map_stencil(m, x, A, b) = [A, " * ", x, " + ", b]
_affine_map_stencil(m, x, A, b::Real) =
    b >= 0 ? [A, " * ", x, " + ", b] :  [A, " * ", x, " - ", abs(b)]

map_stencil_broadcast(m::AbstractAffineMap, x) = _affine_map_stencil_broadcast(m, x, unsafe_matrix(m), unsafe_vector(m))
_affine_map_stencil_broadcast(m, x, A, b) = [A, " .* ", x, " .+ ", b]
_affine_map_stencil_broadcast(m, x, A::Number, b) = [A, " * ", x, " .+ ", b]

map_object_parentheses(m::AbstractAffineMap) = true
map_stencil_parentheses(m::AbstractAffineMap) = true
