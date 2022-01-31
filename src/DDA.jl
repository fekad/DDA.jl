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
abstract type AbstractGrid{T,N} <: AbstractArray{T,N} end
abstract type AbstractShape end
abstract type AbstractPolarizability end
abstract type AbstractScatterer end
abstract type AbstractIncidentField end
abstract type AbstractProblem end
abstract type AbstractMethod end
abstract type AbstractSolution end


export CartesianGrid
include("grid.jl")

export Sphere, Disk
include("shapes.jl")
export CM, CMRR, GOHG, ILDR, LDR
include("polarizability.jl")
include("scatterer.jl")

export PlaneWave, field
include("incidentfield.jl")

export TensorConvolution
include("problems.jl")
include("solvers.jl")

include("farfield.jl")
include("nearfield.jl")

end # module DDA