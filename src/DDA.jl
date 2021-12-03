module DDA

# import Base.CartesianIndecies

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

abstract type AbstractTarget end
# abstract type AbstractTarget{Dim,T<:Number} end
# Base.ndims(::AbstractTarget{Dim}) where {Dim} = Dim


# import Base
export CartesianGrid
include("grid.jl")

export Sphere, Disk
include("targets.jl")


export CM, CMRR, GOHG, ILDR, LDR
include("polarizability.jl")

using Rotations: RotZY
export PlaneWave, field
include("incidentfield.jl")

include("interaction.jl")

include("farfield.jl")

include("nearfield.jl")




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
