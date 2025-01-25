module FunctionMaps

using CompositeTypes, CompositeTypes.Display, CompositeTypes.Indexing
using LinearAlgebra
using StaticArrays

import Base:
    convert, show,
    ==, hash        # for the equality of maps

import CompositeTypes: component, components

# Exhaustive list of exports:

# from util/common.jl
export prectype, numtype

# from generic/map.jl
export Map, MapRef,
    applymap, isequalmap,
    domaintype, codomaintype,
    inverse, leftinverse, rightinverse,
    mapsize, jacobian, jacdet, diffvolume,
    isrealmap
# from generic/canonical.jl
export canonicalmap
# from generic/composite.jl
export composedmap
# from generic/product.jl
export ProductMap, productmap

# from concrete/basic.jl
export IdentityMap,
    ZeroMap,
    UnityMap,
    ConstantMap,
    isconstantmap,
    mapconstant,
    isidentitymap
# from concrete/affine
export AffineMap, Translation, LinearMap,
    affinematrix, affinevector,
    islinearmap, isaffinemap

include("util/common.jl")

include("generic/map.jl")
include("generic/interface.jl")
include("generic/canonical.jl")
include("generic/lazy.jl")
include("generic/inverse.jl")
include("generic/jacobian.jl")
include("generic/composite.jl")
include("generic/product.jl")
include("generic/isomorphism.jl")
include("concrete/basic.jl")
include("concrete/affine/abstractaffine.jl")
include("concrete/affine/linear.jl")
include("concrete/affine/translation.jl")
include("concrete/affine/affine.jl")
include("concrete/affine/special.jl")
include("concrete/coordinates.jl")
include("concrete/arithmetics.jl")

end
