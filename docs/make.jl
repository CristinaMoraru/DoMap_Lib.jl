using DoMap_Lib
using Documenter

DocMeta.setdocmeta!(DoMap_Lib, :DocTestSetup, :(using DoMap_Lib); recursive=true)

makedocs(;
    modules=[DoMap_Lib],
    authors="Cristina Moraru <lilianacristina.moraru@uni-due.de>",
    sitename="DoMap_Lib.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
