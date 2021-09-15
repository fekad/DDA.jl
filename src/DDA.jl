module DDA

using LinearAlgebra
using StaticArrays
using Dierckx
# using Interpolations
# # using Interpolations NOTE: hard to use
# using Core: Vector

import Meshes
import Base.convert
import Base.getindex
import Base.CartesianIndecies

using StaticArrays:SVector


export E_inc, Permittivity, PermittivityTable, permittivity
include("permittivity.jl")
include("polarizability.jl")
include("field.jl")
include("interaction.jl")

export CartesianGrid
include("disrcetisation.jl")

export Sphere, dipoles
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
