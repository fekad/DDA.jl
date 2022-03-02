# 1. load or create the coordinates of the dipoles,
# 2. load or assign the polarizability αj to each dipole,
# 3. calculated the incident field Einc,j at each dipole,
# 4. assemble the interaction matrix A and
# 5. solve for P in the system of linear equations
#

# doi:10.1016/j.jqsrt.2011.03.012
# Knowing Pj, other quantities such as
# - the scattered field,
# - dipole force,
# - Poynting vector,
# - extinction,
# - absorption and
# - scattering cross sections,
# - phase function,
# - Mueller matrix etc.
# can be calculated.



# TODO: allpha is a 3x3 (symmetric) tensor by defualt )
# calc_Ajj(alph) = Diagonal(1 / alph * I, 3)


# """
#     FrequencySimulation([particles::AbstractParticles=[],]
#                         source::AbstractSource)
# Build a FrequencySimulation. If particles are not provided, an empty array is used.
# After building, you can [`run`](@ref) the simulation to get a [`FrequencySimulationResult`](@ref).
# """
# mutable struct FrequencySimulation{Dim,P<:PhysicalMedium} <: Simulation{Dim}
#     "Vector of particles, can be of different types."
#     particles::AbstractParticles
#     "RegularSource wave, where source.medium is the background medium of the simulation."
#     source::AbstractSource{P}
# end

# High level
# mutable struct DipoleProblem <: AbstractProblem
#     "Vector of scatterers, can be of different types."
#     scatterers::Vector{S} where S<:AbstractScatterer
#     "single or as many as scatterers"
#     grid
#     "RegularSource wave, where source.medium is the background medium of the simulation."
#     Einc
# end

# low-level structure (with high level constructors)
struct DipoleProblem <: AbstractProblem
    k
    E0
    dipoles
    alphas
    Eincs
end


# struct DipoleProblem{Dim, T<:AbstractFloat} <: AbstractProblem
#     k::T
#     dipoles::Vector{V} where V<:StaticVector{Dim,T}
#     alphas::Vector{Complex{T}}
#     Eincs::Vector{C} where C<:StaticVector{Dim, <:Complex{T}}
# end


struct DipoleProblemHigh{Dim,T<:AbstractFloat} <: AbstractProblem
    scatterers::Vector{AbstractScatterer} # dipoles + alphas
    Eincs::AbstractIncidentField
end


struct GridProblem <: AbstractProblem
    k
    E0
    grid
    inds
    alphas
    Eincs
end

struct GridProblemFuture <: AbstractProblem
    grid
    scatterers
    Einc
end



# TODO: fix spacing
function polarisbility(m::LDRModel, prob::GridProblem)
    ε = m.ε
    d = prob.grid.spacing[1]
    kvec, E₀  = prob.Einc.kvec, prob.Einc.E₀
    LDR(ε, d, kvec, E₀)
end


