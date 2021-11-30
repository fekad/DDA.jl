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
    origin::SVector{N,T}
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

Base.size(g::CartesianGrid) = g.dims

@inline function Base.getindex(g::CartesianGrid{T,N}, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(g, I...)
    g.origin + (I .- 1) .* g.spacing
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
    print(io, "  spacing: ", Tuple(spacing(g)))
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
