using Test

using LinearAlgebra, StaticArrays
using CompositeTypes, CompositeTypes.Indexing

using FunctionMaps

include("aqua.jl")

include("test_common.jl")
include("test_interface.jl")

include("test_generic.jl")
include("test_basic.jl")
include("test_affine.jl")
include("test_canonical.jl")
include("test_product.jl")
include("test_maps.jl")
include("test_arithmetics.jl")
