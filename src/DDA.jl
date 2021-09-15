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
include("polarizability.jl")
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
