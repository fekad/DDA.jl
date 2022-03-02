# TODO: origin(start) spacing (steps) and dims (length)  VS start stop and length (dims)
# # A regular grid is a tessellation of n-dimensional Euclidean space by congruent parallelotopes (e.g. bricks).
# # A Cartesian grid is a special case where the elements are unit squares or unit cubes, and the vertices are points on the integer lattice.

"""
    CartesianGrid(origin, spacing, dims)

A Cartesian grid with dimensions `dims`, lower left corner at `origin`
and cell spacing `spacing`. The three arguments must have the same length.

    CartesianGrid(dims)
    CartesianGrid(dim1, dim2, ...)

Alternatively, a Cartesian grid can be constructed by only passing the dimensions
`dims` as a tuple, or by passing each dimension `dim1`, `dim2`, ... separately.
In this case, the origin and spacing default to (0,0,...) and (1,1,...).

## Examples

Create a 3D grid with 100x100x50 locations:

```julia
julia> CartesianGrid(100, 100, 50)
```

Create a 2D grid with 100x100 locations and origin at [10.,20.] units:

```julia
julia> CartesianGrid([10.,20.], [1.,1.], (100, 100))
```
"""
struct CartesianGrid{T,N} <: AbstractGrid{SVector{N,T},N}
    origin::Point{N,T}
    spacing::SVector{N,T}
    dims::Dims{N}

    function CartesianGrid{T,N}(origin, spacing, dims) where {T,N}
        @assert all(dims .> 0) "dimensions must be positive"
        @assert all(spacing .> 0) "spacing must be positive"
        new(origin, spacing, dims)
    end
end


# Constructors

function CartesianGrid(origin::AbstractVector{T}, spacing::AbstractVector{T}, dims::Dims{N}) where {T,N}
    CartesianGrid{T,N}(SVector(origin...), SVector(spacing...), dims)
end

function CartesianGrid(origin::AbstractVector, spacing::AbstractVector, dims)
    CartesianGrid(promote(origin, spacing)..., Dims(dims))
end

# Constructors for the unit step grid
function CartesianGrid{T}(dims::Vararg{Int,N}) where {T,N}
    CartesianGrid{T,N}(zero(SVector{N,T}), ones(SVector{N,T}), dims)
end

function CartesianGrid(dims::Vararg{Int,N}) where {N}
    CartesianGrid{Int}(dims...)
end


# AbstractArray interface
# Base.ndims(::CartesianGrid{Dim}) where {Dim} = Dim
Base.size(g::CartesianGrid) = g.dims

@inline function Base.getindex(g::CartesianGrid{T,N}, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(g, I...)
    return @. g.origin + (I - 1) * g.spacing
end


#  Extra functions

Base.minimum(g::CartesianGrid) = g.origin
Base.maximum(g::CartesianGrid) = g.origin + width(g)
Base.extrema(g::CartesianGrid) = (minimum(g), maximum(g))

spacing(g::CartesianGrid) = g.spacing
width(g::CartesianGrid) = (g.dims .- 1) .* g.spacing
center(g::CartesianGrid) = g.origin + width(g) ./ 2

function Base.show(io::IO, g::CartesianGrid{T,N}) where {T,N}
    dims = join(g.dims, "×")
    print(io, "$dims CartesianGrid{$T,$N}")
end

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid)
    println(io, g)
    println(io, "  minimum: ", minimum(g))
    println(io, "  maximum: ", maximum(g))
    print(io, "  spacing: ", spacing(g))
end

# CartesianGrid(origin::SVector{T,N}, spacing::SVector{T,N}, dims::Dims{N}) where {T,N} =
#     CartesianGrid{T,N}(origin, spacing, dims)

# CartesianGrid(origin::NTuple{T,N}, spacing::NTuple{T,N}, dims::Dims{N}) where {T,N} =
#     CartesianGrid{T,N}(SVector(origin), SVector(spacing), dims)

#
# function CartesianGrid{T}(dims::Dims{N}) where {T,N}
#     @show "CartesianGrid{T}(dims::Dims{N}) where {T,N}"
#     CartesianGrid{T,N}(zero(SVector{N,T}), ones(SVector{N,T}), dims)
# end
#
# function CartesianGrid{T}(dims::Vararg{Int,N}) where {T,N}
#     @show "CartesianGrid{T}(dims::Vararg{Int,N}) where {T,N}"
#     CartesianGrid{T}(dims)
# end
#
# function CartesianGrid(dims::Dims{N}) where {N}
#     @show "CartesianGrid(dims::Dims{N}) where {N}"
#     CartesianGrid{Int}(dims)
# end
#
# function CartesianGrid(dims::Vararg{Int,N}) where {N}
#     @show "CartesianGrid(dims::Vararg{Int,N}) where {N}"
#     CartesianGrid{Int}(dims)
# end
#
# function getindex(g::CartesianGrid{T,N}, I::CartesianIndex{N}) where {T,N}
#     coord = g.origin + (I.I .- 1) .* g.spacing
#     return coord
# end

# Future:
# - griddata = domain(grid) + values or inds + values or occ + values
# - grid = topology(dims) + scaling(dx, dy, dz)

#
# # """
# #     CubicGrid(origin, spacing, dims)
# #
# # A Cartesian grid with dimensions `dims`, lower left corner at `origin`
# # and cell spacing `spacing`. The three arguments must have the same length.
# #
# #     CubicGrid(dims)
# #     CubicGrid(dim1, dim2, ...)
# #
# # Alternatively, a Cartesian grid can be constructed by only passing the dimensions
# # `dims` as a tuple, or by passing each dimension `dim1`, `dim2`, ... separately.
# # In this case, the origin and spacing default to (0,0,...) and (1,1,...).
# #
# # ## Examples
# #
# # Create a 3D grid with 100x100x50 locations:
# #
# # ```julia
# # julia> CubicGrid(100, 100, 50)
# # ```
# #
# # Create a 2D grid with 100x100 locations and origin at [10.,20.] units:
# #
# # ```julia
# # julia> CubicGrid([10.,20.], [1.,1.], (100, 100))
# # ```
# # """
# # # AbstractGrid{<:StaticVector{3, T}, 3} where T<:Real
# # # AbstractGrid{V, 3} where {T<:Real,V<:StaticVector{3, T}}
# # # AbstractGrid{V, 3} where V<:StaticVector{3, T} where T<:Real
# # # struct CubicGrida{T<:Real, V<:StaticVector{3,T}} <: AbstractGrid{V,3} where T<:Real
#
# struct CubicGrid{T<:Real, V<:StaticVector{3,T}} <: AbstractGrid{V,3}
#     origin::V
#     spacing::T
#     dims::Dims{3}
#
#     function CubicGrid{T,V}(origin, spacing, dims) where {T<:Real, V<:StaticVector{3,T}}
#         @assert all(dims .> 0) "dimensions must be positive"
#         @assert (spacing > 0) "spacing must be positive"
#         new(origin, spacing, dims)
#     end
# end
# # struct CubicGrid3{T,V,D} where {D,T<:Real,V<:StaticVector{D,T}}
# #     origin::V
# #     spacing::T
# #     dims::Dims{D}
# # end
# #
# # g = CubicGrid(SVector{3}(1,2,3), 1, Dims((10,10,10)))
# # g2 = CubicGrid2(SVector{3}(1,2,3), 1, Dims((10,10,10)))
#
#
#
# #
# #
# # # Constructors
# # function CubicGrid(origin::AbstractVector{V}, spacing::T, dims) where {T<:Real, V<:Real}
# #     TT = promote_type(V,T)
# #     CubicGrid{T,N}(SVector{3,TT}(origin...), TT(spacing), Dims{3}(dims))
# # end
# #
# #
# # # Constructors for the unit step grid
# # function CubicGrid(dims::Vararg{Int,3})
# #     CubicGrid(zero(SVector{3,T}), ones(T), Dims(dims))
# # end
# #
# #
# # # AbstractArray interface
# # Base.size(g::CubicGrid) = g.dims
# #
# # @inline function Base.getindex(g::CubicGrid{T,N}, I::Vararg{Int,N}) where {T,N}
# #     @boundscheck checkbounds(g, I...)
# #     return @. g.origin + (I - 1) * g.spacing
# # end
# #
# #
# # #  Extra functions
# #
# # Base.minimum(g::CubicGrid) = g.origin
# # Base.maximum(g::CubicGrid) = g.origin + width(g)
# # Base.extrema(g::CubicGrid) = (minimum(g), maximum(g))
# #
# # spacing(g::CubicGrid) = g.spacing
# # width(g::CubicGrid) = (g.dims .- 1) .* g.spacing
# # center(g::CubicGrid) = g.origin + width(g) ./ 2
# #
# # function Base.show(io::IO, g::CubicGrid{T,N}) where {T,N}
# #     dims = join(g.dims, "×")
# #     print(io, "$dims CubicGrid{$T,$N}")
# # end
# #
# # function Base.show(io::IO, ::MIME"text/plain", g::CubicGrid)
# #     println(io, g)
# #     println(io, "  minimum: ", minimum(g))
# #     println(io, "  maximum: ", maximum(g))
# #     print(io, "  spacing: ", spacing(g))
# # end
# #
# #