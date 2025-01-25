using FunctionMaps
using Documenter

DocMeta.setdocmeta!(FunctionMaps, :DocTestSetup, :(using FunctionMaps); recursive=true)

makedocs(;
    modules=[FunctionMaps],
    authors="Daan Huybrechs <daan.huybrechs@kuleuven.be> and contributors",
    sitename="FunctionMaps.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
	    "Examples" => "examples.md",
        "API" => Any[
            "Public API Reference" => "api.md",
            "Internal API Reference" => "internal.md"
        ],
    ],
)

deploydocs(;
    repo="github.com/JuliaApproximation/FunctionMaps.jl.git", devbranch="main",
)
