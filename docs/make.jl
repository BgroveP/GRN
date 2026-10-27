using Documenter
using GRN

makedocs(
    sitename = "GRN",
    modules = [GRN],
    format = Documenter.HTML(),
    pages = ["Home" => "index.md"],
    repo = "github.com/BgroveP/GRN/blob/{commit}{path}#{line}",
)

deploydocs(repo = "github.com/BgroveP/GRN.git", push_preview = true)
