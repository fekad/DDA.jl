module DDA

using LinearAlgebra
using StaticArrays
using Dierckx
# using Interpolations
# # using Interpolations NOTE: hard to use
# using Core: Vector


export E_inc, Permittivity, PermittivityTable, permittivity
include("permittivity.jl")
include("polarizability.jl")
include("fields.jl")
include("dipole.jl")
include("disrcetisation.jl")

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
