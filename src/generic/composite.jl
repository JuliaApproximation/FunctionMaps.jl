# ensure that the domaintype of each map in the sequence agrees with the codomain type of the
# preceding map
function match_domain_codomain_types(::Type{T}, map, maps...) where {T}
    Y = promote_type(T, domaintype(map))
    m1 = convert_domaintype(Y, map)
    (m1, match_domain_codomain_types(codomaintype(m1), maps...)...)
end

match_domain_codomain_types(::Type{T}, map) where {T} = (convert_domaintype(T, map),)

"""
    ComposedMap{T,MAPS}

The composition of several maps.

The `components` of a `ComposedMap` are the maps in the order that they are applied
to the input.
"""
struct ComposedMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

function ComposedMap(map1, maps...)
    P = prectype(map1, maps...)
    if P == Any
        # don't try to promote types
        _ComposedMap(domaintype(map1), map1, maps...)
    else
        T = to_prectype(P, domaintype(map1))
        _ComposedMap(T, match_domain_codomain_types(T, map1, maps...)...)
    end
end
ComposedMap{T}(map1, maps...) where {T} = _ComposedMap(T, match_domain_codomain_types(T, map1, maps...)...)
_ComposedMap(::Type{T}, maps...) where {T} = ComposedMap{T,typeof(maps)}(maps)

similarmap(m::ComposedMap, ::Type{T}) where {T} = ComposedMap{T}(components(m)...)

codomaintype(m::ComposedMap) = codomaintype(last(m.maps))

# Maps are applied in the order that they appear in m.maps
applymap(m::ComposedMap, x) = applymap_rec(x, m.maps...)
applymap_rec(x) = x
applymap_rec(x, map1, maps...) = applymap_rec(applymap(map1, x), maps...)

# The size of a composite map depends on the first and the last map to be applied
# We check whether they are scalar_to_vector, vector_to_vector, etcetera
mapsize(m::ComposedMap) = _composed_mapsize(m, last(m.maps), first(m.maps), mapsize(last(m.maps)), mapsize(first(m.maps)))
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int,Int}, S1::Tuple{Int,Int}) = (S_end[1],S1[2])
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int,Int}, S1::Tuple{Int}) =
    is_vector_to_scalar(m_end) ? () : (S_end[1],)
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int,Int}, S1::Tuple{}) =
    is_vector_to_scalar(m_end) ? () : (S_end[1],)
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int}, S1::Tuple{Int,Int}) = (S_end[1],S1[2])
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int}, S1::Tuple{Int}) = (S_end[1],)
_composed_mapsize(m, m_end, m1, S_end::Tuple{Int}, S1::Tuple{}) = (S_end[1],)
_composed_mapsize(m, m_end, m1, S_end::Tuple{}, S1::Tuple{Int,Int}) = (1,S1[2])
_composed_mapsize(m, m_end, m1, S_end::Tuple{}, S1::Tuple{Int}) = ()
_composed_mapsize(m, m_end, m1, S_end::Tuple{}, S1::Tuple{}) = ()

function jacobian(m::ComposedMap, x)
    f, fd = backpropagate(x, reverse(components(m))...)
    fd
end
backpropagate(x, m1) = (applymap(m1, x), jacobian(m1, x))
function backpropagate(x, m2, ms...)
    f, fd = backpropagate(x, ms...)
    applymap(m2, f), jacobian(m2, f) * fd
end

for op in (:inverse, :leftinverse, :rightinverse)
    @eval $op(cmap::ComposedMap) = ComposedMap(reverse(map($op, components(cmap)))...)
end

inverse(m::ComposedMap, x) = inverse_rec(x, reverse(components(m))...)
inverse_rec(x) = x
inverse_rec(x, map1, maps...) = inverse_rec(inverse(map1, x), maps...)

leftinverse(m::ComposedMap, x) = leftinverse_rec(x, reverse(components(m))...)
leftinverse_rec(x) = x
leftinverse_rec(x, map1, maps...) = leftinverse_rec(leftinverse(map1, x), maps...)

rightinverse(m::ComposedMap, x) = rightinverse_rec(x, reverse(components(m))...)
rightinverse_rec(x) = x
rightinverse_rec(x, map1, maps...) = rightinverse_rec(rightinverse(map1, x), maps...)


promote_composing_maps(m1, m2) =
    _promote_composing_maps(codomaintype(m1), domaintype(m2), m1, m2)
function _promote_composing_maps(::Type{S}, ::Type{T}, m1, m2) where {S,T}
    U = promote_type(S, T)
    convert_codomaintype(U, m1), convert_domaintype(U, m2)
end

composedmap() = ()
composedmap(m) = m
composedmap(m1, m2) = composedmap1(promote_composing_maps(m1, m2)...)
composedmap1(m1, m2) = composedmap2(m1, m2)
composedmap2(m1, m2) = default_composedmap(m1, m2)
default_composedmap(m1, m2) = ComposedMap(m1, m2)

composedmap(m1, m2, maps...) = composedmap(composedmap(m1, m2), maps...)

composedmap(m1::ComposedMap, m2::ComposedMap) =
    ComposedMap(components(m1)..., components(m2)...)
function composedmap1(m1::ComposedMap, m2)
    T = promote_type(codomaintype(m1), domaintype(m2))
    ComposedMap(components(m1)..., convert_domaintype(T, m2))
end
function composedmap2(m1, m2::ComposedMap)
    T = promote_type(codomaintype(m1), domaintype(m2))
    ComposedMap(convert_codomaintype(T, m1), components(m2)...)
end

# Arguments to ∘ should be reversed before passing on to mapcompose
Base.:∘(map1::Map, map2::Map) = composedmap(map2, map1)


isequalmap(m1::ComposedMap, m2::ComposedMap) =
    ncomponents(m1) == ncomponents(m2) && all(map(isequalmap, components(m1), components(m2)))
map_hash(m::ComposedMap, h::UInt) = hashrec("ComposedMap", collect(components(m)), h)

Display.combinationsymbol(m::ComposedMap) = Display.Symbol('∘')
Display.displaystencil(m::ComposedMap) =
    composite_displaystencil(m; reversecomponents=true)
map_object_parentheses(m::ComposedMap) = true
map_stencil_parentheses(m::ComposedMap) = true
show(io::IO, mime::MIME"text/plain", m::ComposedMap) = composite_show(io, mime, m)
show(io::IO, m::ComposedMap) = composite_show_compact(io, m)

## Lazy multiplication

"The lazy multiplication of one or more maps."
struct MulMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

MulMap(map1::Map, maps::Map...) = MulMap(promote_maps(map1, maps...)...)
MulMap(map1::Map{T}, maps::Map{T}...) where {T} = MulMap{T}(map1, maps...)
MulMap{T}(maps::Map{T}...) where {T} = MulMap{T,typeof(maps)}(maps)
MulMap{T}(maps...) where {T} = _mulmap(T, convert_domaintype.(Ref(T), maps)...)
_mulmap(::Type{T}, maps...) where {T} = MulMap{T,typeof(maps)}(maps)

similarmap(m::MulMap, ::Type{T}) where {T} = MulMap{T}(m.maps...)

applymap(m::MulMap, x) = reduce(*, applymap.(components(m), Ref(x)))

multiply_map() = ()
multiply_map(m) = m
multiply_map(m1, m2) = multiply_map1(m1, m2)
multiply_map1(m1, m2) = multiply_map2(m1, m2)
multiply_map2(m1, m2) = default_multiply_map(m1, m2)
default_multiply_map(m1, m2) = MulMap(m1, m2)

multiply_map(m1, m2, maps...) = multiply_map(multiply_map(m1, m2), maps...)

multiply_map(m1::MulMap, m2::MulMap) =
    MulMap(components(m1)..., components(m2)...)
multiply_map1(m1::MulMap, m2) = MulMap(components(m1)..., m2)
multiply_map2(m1, m2::MulMap) = MulMap(m1, components(m2)...)

function Display.displaystencil(m::MulMap)
    A = Any[]
    list = components(m)
    push!(A, "x -> ")
    push!(A, Display.SymbolObject(list[1]))
    push!(A, "(x)")
    for i in 2:length(list)
        push!(A, " * ")
        push!(A, Display.SymbolObject(list[i]))
        push!(A, "(x)")
    end
    A
end
show(io::IO, mime::MIME"text/plain", m::MulMap) = composite_show(io, mime, m)


## Lazy sum

"The lazy sum of one or more maps."
struct SumMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

SumMap(map1::Map, maps::Map...) = SumMap(promote_maps(map1, maps...)...)
SumMap(map1::Map{T}, maps::Map{T}...) where {T} = SumMap{T}(map1, maps...)
SumMap{T}(maps::Map{T}...) where {T} = SumMap{T,typeof(maps)}(maps)
SumMap{T}(maps...) where {T} = _summap(T, convert_domaintype.(Ref(T), maps)...)
_summap(::Type{T}, maps...) where {T} = SumMap{T,typeof(maps)}(maps)

similarmap(m::SumMap, ::Type{T}) where {T} = SumMap{T}(m.maps...)

applymap(m::SumMap, x) = reduce(+, applymap.(components(m), Ref(x)))

sum_map() = ()
sum_map(m) = m
sum_map(m1, m2) = sum_map1(m1, m2)
sum_map1(m1, m2) = sum_map2(m1, m2)
sum_map2(m1, m2) = default_sum_map(m1, m2)
default_sum_map(m1, m2) = SumMap(m1, m2)

sum_map(m1, m2, maps...) = sum_map(sum_map(m1, m2), maps...)

sum_map(m1::SumMap, m2::SumMap) =
    SumMap(components(m1)..., components(m2)...)
sum_map1(m1::SumMap, m2) = SumMap(components(m1)..., m2)
sum_map2(m1, m2::SumMap) = SumMap(m1, components(m2)...)

function Display.displaystencil(m::SumMap)
    A = Any[]
    list = components(m)
    push!(A, "x -> ")
    push!(A, Display.SymbolObject(list[1]))
    push!(A, "(x)")
    for i in 2:length(list)
        push!(A, " + ")
        push!(A, Display.SymbolObject(list[i]))
        push!(A, "(x)")
    end
    A
end
show(io::IO, mime::MIME"text/plain", m::SumMap) = composite_show(io, mime, m)


# Define the jacobian of a composite map
jacobian(m::ComposedMap) = composite_jacobian(reverse(components(m))...)
composite_jacobian(map1) = jacobian(map1)
composite_jacobian(map1, map2) = multiply_map(composedmap(map2, jacobian(map1)), jacobian(map2))
function composite_jacobian(map1, map2, maps...)
    rest = ComposedMap(reverse(maps)..., map2)
    f1 = composedmap(rest, jacobian(map1))
    f2 = composite_jacobian(map2, maps...)
    multiply_map(f1, f2)
end

jacobian(m::MulMap) = mul_jacobian(components(m)...)
mul_jacobian() = ()
mul_jacobian(map1) = jacobian(map1)
mul_jacobian(map1, map2) = sum_map(multiply_map(jacobian(map1), map2), multiply_map(map1, jacobian(map2)))
function mul_jacobian(map1, map2, maps...)
    rest = multiply_map(map2, maps...)
    mul_jacobian(map1, rest)
end
function jacobian(m::MulMap, x)
    z = map(t -> applymap(t,x), components(m))
    zd = map(t -> jacobian(t, x), components(m))
    sum(prod(z[1:i-1]) * zd[i] * prod(z[i+1:end]) for i in 1:ncomponents(m))
end


jacobian(m::SumMap) = sum_jacobian(components(m)...)
sum_jacobian() = ()
sum_jacobian(map1) = jacobian(map1)
sum_jacobian(maps...) = sum_map(map(jacobian, maps)...)

jacobian(m::SumMap, x) = sum(jacobian(mc, x) for mc in components(m))
