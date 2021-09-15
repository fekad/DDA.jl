# DDSCAT approach:
# - lattice (basis) vector
# - lattice spacings
# - lattice offset
# - integer indecies
#
# DDSCAT example
#  >TARELL  ellipsoidal grain; AX,AY,AZ= 30.0000 30.0000 30.0000
# 14328 = NAT
# 1.000000  0.000000  0.000000 = A_1 vector
# 0.000000  1.000000  0.000000 = A_2 vector
# 1.000000  1.000000  1.000000 = lattice spacings (d_x,d_y,d_z)/d
#  0.50000   0.50000   0.50000 = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d for dipole 0 0 0
#    JA  IX  IY  IZ ICOMP(x,y,z)
#     1   -2   -4  -15 1 1 1

# using Meshes: CartesianGrid
# This implemetation is based on the CartesianGrid in Meshes.jl
# Diferences:
# - the elements are the centroids of the original grid elelments
# - iterataion over coordinates (using CartesianIndecies)
# - simplified topology

abstract type AbstractGrid{Dim,T} end

struct CartesianGrid{Dim,T} <: AbstractGrid{Dim,T}
    origin::SVector{Dim,T}
    spacing::SVector{Dim,T}
    dims::Dims{Dim}

    function CartesianGrid{Dim,T}(origin, spacing, dims) where {Dim,T}
        @assert all(dims .> 0) "dimensions must be positive"
        @assert all(spacing .> 0) "spacing must be positive"
        new(origin, spacing, dims)
    end
end

# Constructors

CartesianGrid(origin::SVector{Dim,T}, spacing::SVector{Dim,T}, dims::Dims{Dim}) where {Dim,T} =
    CartesianGrid{Dim,T}(origin, spacing, dims)

CartesianGrid(origin::NTuple{Dim,T}, spacing::NTuple{Dim,T}, dims::Dims{Dim}) where {Dim,T} =
    CartesianGrid{Dim,T}(SVector(origin), SVector(spacing), dims)

CartesianGrid(origin::AbstractVector{T}, spacing::AbstractVector{T}, dims::Dims{Dim}) where {Dim,T} =
    CartesianGrid(SVector(origin...), SVector(spacing...), Dims(dims))

CartesianGrid(origin::AbstractVector, spacing::AbstractVector, dims::Dims{Dim}) where {Dim} =
    CartesianGrid(promote(origin, spacing)..., Dims(dims))

# Constructors for unit step grid

CartesianGrid{T}(dims::Dims{Dim}) where {Dim,T} =
CartesianGrid{Dim,T}(zero(SVector{Dim,T}), ones(SVector{Dim,T}), dims)
CartesianGrid{T}(dims::Vararg{Int,Dim}) where {Dim,T} = CartesianGrid{T}(dims)

CartesianGrid(dims::Dims{Dim}) where {Dim} = CartesianGrid{Int}(dims)
CartesianGrid(dims::Vararg{Int,Dim}) where {Dim} = CartesianGrid{Int}(dims)


# Indexing functions

Base.CartesianIndices(g::CartesianGrid{Dim}) where {Dim} = CartesianIndices(g.dims)
function Base.getindex(g::CartesianGrid{Dim}, I::CartesianIndex{Dim}) where {Dim}
    coord  = g.origin + (I.I .- 1) .* g.spacing
    return coord
end

size(g::CartesianGrid) = g.dims

# Simplified version of CartesianGrid
struct CubicGrid
    origin::Vector{Float64}
    spacing::Vector{Float64}
    dims::Dims{3}

    function CubicGrid(origin, spacing, dims)
        @assert all(dims .> 0) "dimensions must be positive"
        @assert all(spacing .> 0) "spacing must be positive"
        # @show origin, spacing, dims
        new(origin, spacing, dims)
    end
end

CubicGrid(dims::Vararg{Int,3}) = CubicGrid(zeros(3), ones(3), Dims(dims))

Meshes.CartesianGrid(g::CubicGrid) = Meshes.CartesianGrid{3,Float64}(g.dims, g.origin .- g.spacing / 2, g.spacing)
convert(::Type{Meshes.CartesianGrid}, g::CubicGrid) = Meshes.CartesianGrid(g)


@inline Base.@propagate_inbounds function Base.getindex(
    g::CubicGrid,
    i::Vararg{Int,3},
)
    #  @boundscheck checkbounds(g, i...)
    return @. g.origin + (i - 1) * g.spacing
end

# Alternative implementations:

# # more abstract
# struct LinRange{T,L<:Integer} <: AbstractRange{T}
#     start::T
#     stop::T
#     len::L
# end

# # using Base: Float64, disable_library_threading_hooks
# using LazyGrids:ndgrid
#
# abstract type AbstractGrid end
#
# import Base:length
#
# struct CubicGrid <: AbstractGrid
#     xrange::AbstractRange
#     yrange::AbstractRange
#     zrange::AbstractRange
#     function CubicGrid(xrange, yrange, zrange)
#         @assert isapprox(step(xrange), step(yrange))
#         @assert isapprox(step(xrange), step(zrange))
#         return new(xrange, yrange, zrange)
#     end
# end

# CubicGrid(xmin, xmax, ymin, ymax, zmin, zmax, dx) = CubicGrid(range(xmin, xmax, step=dx), range(ymin, ymax, step=dx), range(zmin, zmax, step=dx))
# CubicGrid(xmin, xmax, ymin, ymax, zmin, zmax, d) = ndgrid(xmin:d:xmax, ymin:d:ymax, zmin:d:zmax)

# dims(g::CubicGrid) = (length(g.xrange), length(g.yrange), length(g.zrange))
# length(g::CubicGrid) = *(dims(g)...)


# cubic grid
# - step size (dx == dy == dz)
# - start(3) stop(3) step()
# - start(3) stop(3)

# struct CubicGrid <: AbstractGrid
#     xi::StepRangeLen
#     yi::StepRangeLen
#     zi::StepRangeLen
# end
# origin(g::CubicGrid) = [g.xi[1], g.yi[1], g.zi[1]]
# dims(g::CubicGrid) = (length(g.xi), length(g.yi), length(g.zi))

# struct CubicLattice <: AbstractGrid
#     origin::SVector{3,Float64}
#     dims::SVector{3,Int64}
#     spacing::Float64
# end
#
# struct RotatadCubicLattice <: AbstractGrid
#     start::SVector{3,Float64}
#     length::SVector{3,Int64}
#     step::Float64
#     orientation::SMatrix{3,3,Float64}
# end
#

# # https://juliageometry.github.io/Meshes.jl/stable/meshes.html#Meshes.CartesianGrid
# # https://github.com/JuliaGeometry/Meshes.jl/blob/master/src/mesh/cartesiangrid.jl
# # https://juliaarrays.github.io/LazyGrids.jl/stable/examples/1-ndgrid/#D-case
#
#
# module StructuredGrids
#
# import Base: axes, getindex, size
#
# using Base: tail, @propagate_inbounds
#
# export grid
#
# struct Grid{T,N,RT} <: AbstractArray{T,N}
#     ranges::RT
# end
#
# @inline eltypes(::Tuple{}) = ()
# @inline eltypes(ranges) = (eltype(first(ranges)), eltypes(tail(ranges))...)
#
# @inline _axes(::Tuple{}) = ()
# @inline _axes(ranges) = (axes(first(ranges), 1), _axes(tail(ranges))...)
# @inline axes(g::Grid) = _axes(g.ranges)
#
# @inline _size(::Tuple{}) = ()
# @inline _size(ranges) = (size(first(ranges), 1), _size(tail(ranges))...)
# @inline size(g::Grid) = _size(g.ranges)
#
#
# _getindex(::Tuple{}, ::Any) = ()
# @propagate_inbounds _getindex(ranges, I) = (first(ranges)[first(I)], _getindex(tail(ranges), tail(I))...)
# @propagate_inbounds getindex(g::Grid, I::Vararg{Int,N}) where {N} = _getindex(g.ranges, I)
#
# grid(ranges::Vararg{Any,N}) where {N} = Grid{Tuple{eltypes(ranges)...},N,typeof(ranges)}(ranges)
#
# end # module