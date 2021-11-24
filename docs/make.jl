using Documenter
using DDA

makedocs(
    sitename = "DDA",
    format = Documenter.HTML(),
    modules = [DDA]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
