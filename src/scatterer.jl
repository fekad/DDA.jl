abstract type AbstractScatterer end

# TODO: Handle multiple targets

struct Scatterer{T<:AbstractShape,M<:AbstractPolarizability} <: AbstractScatterer
    target::T
    model::M
end
