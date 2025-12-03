using RayTracingEngine
using Documenter

DocMeta.setdocmeta!(RayTracingEngine, :DocTestSetup, :(using RayTracingEngine); recursive=true)

makedocs(;
    modules=[RayTracingEngine],
    authors="F. Capolupo",
    sitename="RayTracingEngine.jl",
    format=Documenter.HTML(; canonical="https://FraCpl.github.io/RayTracingEngine.jl", edit_link="master", assets=String[]),
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/FraCpl/RayTracingEngine.jl", devbranch="master")
