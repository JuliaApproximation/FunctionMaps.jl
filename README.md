# FunctionMaps

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaApproximation.github.io/FunctionMaps.jl/dev)
[![Build Status](https://github.com/JuliaApproximation/FunctionMaps.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaApproximation/FunctionMaps.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaApproximation/FunctionMaps.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/FunctionMaps.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

FunctionMaps.jl is a package designed to represent scalar- and vector-valued functions. It makes it easy to make new maps from existing ones, via composition or by taking products.

## Examples

For more information, see the [documentation](https://JuliaApproximation.github.io/FunctionMaps.jl/dev).

### Linear and affine maps

In comparison to existing packages with similar functionality, FunctionMaps aims to support non-invertible maps. These can be used to map between spaces of different dimension, such as from 2D to 3D. This leads to a linear map involving a rectangular matrix:
```julia
julia> using FunctionMaps

julia> m = LinearMap([0 0; 1 0; 0 1.0])
x -> A * x

A = 3×2 Matrix{Float64}:
 0.0  0.0
 1.0  0.0
 0.0  1.0

julia> m([1; 2])
3-element Vector{Float64}:
 0.0
 1.0
 2.0
```
This function `m` maps the 2D plane to the YZ-plane in a three-dimensional Euclidean space.

The map is not invertible on the full 3D space, but it is invertible on its range. In other words, it does have a left inverse from the YZ-plane back to 2D, which is again represented by a rectangular matrix:
```julia
julia> mi = leftinverse(m)
x -> A * x

A = 2×3 Matrix{Float64}:
 0.0  1.0  0.0
 0.0  0.0  1.0

julia> mi(m([1; 2]))
2-element Vector{Float64}:
 1.0
 2.0
```
The behaviour of the left inverse outside the image of the map is not defined in general, as indeed there may be multiple left inverses satisfying `mi(m(x)) == x`, but in this case it corresponds to a projection onto the YZ-plane.
```julia
julia> mi([1; 2; 3])
2-element Vector{Float64}:
 2.0
 3.0
```

### Domain and codomain type

Maps tend to be flexible in the type of the argument if the call makes sense mathematically. Yet, most maps do support one particular domain type. In that sense, in Julia terminology maps more closely resemble a method than a generic function.

Maps are functions of a single variable only. Functions of several variables can be represented as vector-valued maps. Such is the case in the current example.
```julia
julia> domaintype(m)
Vector{Float64} (alias for Array{Float64, 1})

julia> codomaintype(m)
Vector{Float64} (alias for Array{Float64, 1})
```


## Maps as an interface

The concrete maps defined in this package inherit from the abstract supertype `Map{T}`. Here, `T` is the domaintype and any instance conceptually represents a function `f(x::T)`. Like Julia functions and methods, the output is not explicitly typed.

A lot of the functionality in the package does not require maps to inherit from `Map{T}`, as long as they behave like maps. They should implement the following interface.

| Required methods | Brief description | Default behaviour |
| ---------------- | ----------------- | ----------------- |
| `applymap(m, x)` | Applies the map `m` to the argument `x` | `m(x)` |

Optional methods include:

| Important optional methods | Default definition | Brief description
| --- | --- | --- |
| `domaintype(m)` | `Any` | Returns a valid or intended type for arguments of the map |
| `codomaintype(m[, ::Type{T}])` | `Any` | The return type of the map (for arguments of type `T`) |


## History and development

The functionality in this package was originally developed as part of the [DomainSets.jl](https://github.com/JuliaApproximation/DomainSets.jl) package. It was extracted to allow more focused development and a better integration with existing packages with similar functionality.

