module DDA

# import Base.CartesianIndecies

using FFTW
using AbstractFFTs
using IterativeSolvers

import Base: *

# Strandrd libraries
using LinearAlgebra
using AbstractFFTs

using StaticArrays
# using StaticArrays:SVector

using GeometryBasics: Point, Point3
export Point, Point3

abstract type AbstractGrid{T,N} <: AbstractArray{T,N} end
abstract type AbstractField end
abstract type AbstractPoint{Dim,T} <: StaticVector{Dim,T} end


# import Base
export CartesianGrid
include("grid.jl")

export Sphere, Disk
include("shapes.jl")


export CM, CMRR, GOHG, ILDR, LDR
include("polarizability.jl")

using Rotations: RotZY
export PlaneWave, field
include("incidentfield.jl")

include("interaction.jl")

include("farfield.jl")

include("nearfield.jl")

export TensorConvolution
include("convolution.jl")

# interfaces for the solver
abstract type PolirazabilityModels end

struct LDRModel <: PolirazabilityModels
    ε
end

# TODO: Handle multiple targets
struct Scatterer
    target
    model
end



abstract type AbstractMethod end

struct Direct <: AbstractMethod end
struct BiCGStabl <: AbstractMethod end
struct BiCGStablFFT <: AbstractMethod end

# CG
# MINRES
# GMRES
# IDRs

abstract type AbstractProblem end

# todo: a_eff
struct GridProblem <: AbstractProblem
    grid
    scatterer
    E_inc
end

struct DipoleProblem <: AbstractProblem
    dipoles
    a_eff
    scatterer
    E_inc
end

# TODO: Implement the following !!!
# α = LDR(ε, d, kvec, E₀)
# alphas = fill(α, length(r2))
polarisbility(target, prob) = 1
polarisbility(target, prob) = 1

# TODO: make it unitless (x=k*a) (Works only with spheres)
function solve(p::GridProblem, alg::BiCGStablFFT;
    reltol=1e-3, verbose=true, kwargs...)
    # @show p kwargs

    # 1. create the coordinates of the dipoles,
    coords, occ = DDA.dipoles(p.grid, p.scatterer.target)

    # 2. assign the polarizability αj to each dipole
    alphas = polarisbility(p.scatterer.model, p)

    # all in one go
    # coords, alphas, occ = DDA.discretize(p.scatterer, p.grid, p.E_inc)



    # 3. calculated the incident field E_inc, at each dipole,
    E_inc = field(p.E_inc, coords)

    k = norm(p.E_inc.k)

    # 4. assemble the interaction matrix A and

    A_conv = TensorConvolution(g, occ, k, alphas)

    # 5. solve for P in the system of linear equations

    Ei = reinterpret(ComplexF64, E_inc)
    P = bicgstabl(A_conv, Ei; reltol=reltol, verbose=verbose)

    P = reinterpret(SVector{3,ComplexF64}, P3)



    return GridSolution(P, prob, alg)
end


abstract type AbstractSolution end
# TODO: storing the alg and its parameters
# TODO: storing the prepared (factorised) interaction matrix
# TODO: do not trim A*E because it is useful for the nearfield calculations
#
struct GridSolution <: AbstractSolution
    P
    prob
    alg
end





# Dipole coordinates
# Polarizability
# Solver
# - incident field
# - interaction matrix
# - extinction scattering and absorption cress sections
# - scattered and total field
#


# Challanges:
# Symmetric{BlockMatrix{Symmetric{Matrix{Int64}}}}
# a = Symmetric(BlockArray(undef_blocks, Symmetric{Int64, Matrix{Int64}}, [2,2], [2,2]))


end # module
