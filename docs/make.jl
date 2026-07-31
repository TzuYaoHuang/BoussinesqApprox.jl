using BoussinesqApprox
using Documenter

DocMeta.setdocmeta!(BoussinesqApprox, :DocTestSetup, :(using BoussinesqApprox); recursive=true)

makedocs(;
    modules=[BoussinesqApprox],
    authors="Tzu-Yao Huang <tzuyao.jason.huang@gmail.com>",
    sitename="BoussinesqApprox.jl",
    format=Documenter.HTML(;
        canonical="https://TzuYaoHuang.github.io/BoussinesqApprox.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/TzuYaoHuang/BoussinesqApprox.jl",
    devbranch="main",
)
