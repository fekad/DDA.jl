module DDA

import Base: size, getindex, convert
# import Base.CartesianIndecies

# Strandrd libraries
using LinearAlgebra
using AbstractFFTs

using StaticArrays
# using StaticArrays:SVector

abstract type AbstractGrid{Dim,T} end
abstract type AbstractField end
abstract type AbstractTarget end


export CartesianGrid
include("disrcetisation.jl")

export Sphere, Disk
include("targets.jl")


# export E_inc, Permittivity, PermittivityTable, permittivity
include("MaterialModels.jl")





export CM, LDR
include("polarizability.jl")

using Rotations: RotZY
export PlaneWave, field
include("incidentfield.jl")

include("interaction.jl")




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
