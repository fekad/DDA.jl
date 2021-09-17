module DDA


import Base: size, getindex, convert
# import Base.CartesianIndecies

# Strandrd libraries
using StaticArrays
using StaticArrays:SVector

using LinearAlgebra

# Interpolation
using Dierckx
# using Interpolations

import Meshes


export E_inc, Permittivity, PermittivityTable, permittivity
include("permittivity.jl")

export polarizability_CM, polarizability_LDR
include("polarizability.jl")

export PlaneWave, field
include("field.jl")

include("interaction.jl")

export CartesianGrid
include("disrcetisation.jl")

export Disk, Sphere, dipoles, discretize
include("scatterer.jl")

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
