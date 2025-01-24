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
    ],
)
