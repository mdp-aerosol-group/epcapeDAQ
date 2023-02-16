using Documenter

makedocs(
  sitename = "EPCAPE-CCN",
  authors = "Markus Petters",
  pages = Any[
    "Home" => "index.md",
    "Quicklooks" => "quicklook.md",
    "Housekeeping" => "housekeeping.md",
  ]
)