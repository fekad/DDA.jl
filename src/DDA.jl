"""
    DDA.jl - Discrete Dipole Approximation

"""
module DDA

using Reexport

using StaticArrays
using LinearAlgebra
import Base: *

using IterativeSolvers
using AbstractFFTs
using FFTW

using GeometryBasics: Point, Point3
export Point, Point3

using Rotations: RotZY
using SpecialFunctions: besselk
using PhysicalConstants.CODATA2018: ħ, ε_0, e as q

# Definition of interfaces
# abstract type AbstractGrid{T,N}end
abstract type AbstractGrid{T,N} <: AbstractArray{T,N} end
abstract type AbstractShape end
abstract type AbstractPolarizability end
abstract type AbstractScatterer end
abstract type AbstractIncidentField end
abstract type AbstractProblem end
abstract type AbstractMethod end
abstract type AbstractSolution end


export CartesianGrid
# export spacing, width, center
include("grid.jl")

export Dipole, Sphere, Disk
export discretize
include("shapes.jl")

export  Atomic, LDRModel
export CM, CMRR, GOHG, ILDR, LDR
include("polarizability.jl")
export Scatterer
include("scatterer.jl")

export PlaneWave, field
include("incidentfield.jl")

export C_abs, C_ext, C_sca
include("crossection.jl")
include("electricfield.jl")

export GridProblem, DipoleProblem
include("problems.jl")

export solve, Direct, BiCGStabl, BiCGStablFFT
include("solvers.jl")


end # module DDA