using Pkg; Pkg.activate("docs")
using Documenter
using AmphiDEB

makedocs(
    sitename = "AmphiDEB", 
    format = Documenter.HTML()
    )

deploydocs(
    repo = "github.com/SimonHansul/AmphiDEB.git",
)