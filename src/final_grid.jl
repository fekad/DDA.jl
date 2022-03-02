using StaticArrays

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
struct CartesianGrid1{Dim,T<:AbstractFloat,V<:StaticVector{Dim,T}} <: AbstractArray{V,Dim}
    origin::V
    spacing::T
    dims::Dims{Dim}

    function CartesianGrid1{Dim,T,V}(origin, spacing, dims) where {Dim,T,V}
        @assert all(dims .> 0) "dimensions must be positive"
        @assert spacing > 0 "spacing must be positive"
        new(origin, spacing, dims)
    end
end

function CartesianGrid1(origin::StaticVector{N,T1}, spacing::T2, dims::Dims{N}) where {T1<:AbstractFloat,T2<:AbstractFloat,N}
    T = promote_type(T1, T2)
    CartesianGrid1(convert.(T, origin), convert(T, spacing), dims)
end

# function CartesianGrid1(origin::V, spacing::T2, dims::Dims{N}) where {N, T1<:AbstractFloat, V<:StaticVector{N,T1},T2<:AbstractFloat}
#     T = promote_type(T1, T2)
#     CartesianGrid1{N,T,V}(convert.(T, origin), convert(T, spacing), dims)
# end

# AbstractVector but not StaticVector -> use Svector as fallback for default value
function CartesianGrid1(origin::AbstractVector{T1}, spacing::T2, dims::Dims{N}) where {T1<:AbstractFloat,T2<:AbstractFloat,N}
    @assert length(origin) == N "Length of origin must match number of dimensions"
    T = promote_type(T1, T2)
    CartesianGrid1(SVector{N,T}(origin...), convert(T, spacing), dims)
end

Base.size(g::CartesianGrid1) = g.dims
Base.ndims(::CartesianGrid1{Dim}) where {Dim} = Dim

Base.minimum(g::CartesianGrid1) = g.origin
Base.maximum(g::CartesianGrid1) = g.origin + width(g)
Base.extrema(g::CartesianGrid1) = (minimum(g), maximum(g))

spacing(g::CartesianGrid1) = g.spacing
width(g::CartesianGrid1{Dim,T,V}) where {Dim,T,V} = V((g.dims .- 1) .* g.spacing)
center(g::CartesianGrid1) = g.origin + width(g) / 2

@inline function Base.getindex(g::CartesianGrid1{Dim}, I::Vararg{Int,Dim}) where {Dim}
    @boundscheck checkbounds(g, I...)
    return @. g.origin + (I - 1) * g.spacing
end

function Base.show(io::IO, g::CartesianGrid1{Dim,T,V}) where {Dim,T,V}
    dims = join(g.dims, "×")
    print(io, "$dims CartesianGrid1{$Dim,$T,$V}")
end

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid1)
    println(io, g)
    println(io, "  minimum: ", minimum(g))
    println(io, "  maximum: ", maximum(g))
    print(io, "  spacing: ", spacing(g))
end



g1 = CartesianGrid1(SVector{3,Float64}(1, 2, 3), 1., Dims((10, 10, 10)))
g2 = CartesianGrid1(SVector(1, 2, 3.0), 1.0, (10, 10, 10))
g3 = CartesianGrid1(SVector(1, 2, 3), 1.0, (10, 10, 10))  # missmatch type (promotion)
g4 = CartesianGrid1([1, 2, 3], 1.0, (10, 10, 10)) # Svector fallback as default
CartesianGrid1([1, 2, 3], 1.0, (3, 3, 3))




# # A regular grid is a tessellation of n-dimensional Euclidean space by congruent parallelotopes (e.g. bricks).
# # A Cartesian grid is a special case where the elements are unit squares or unit cubes, and the vertices are points on the integer lattice.

# TODO: origin(start) spacing (steps) and dims (length)  VS start stop and length (dims)
# abstract type AbstractGrid{Dim,V} end
# abstract type StructuredGrid{Dim,V} <: AbstractGrid{Dim,V} end
# abstract type UnstructuredGrid{Dim,V} <: AbstractGrid{Dim,V} end

# CartesianGrid{T,D} <: StructuredGrid{T,D} (spacing, dims)
# RegularGrid{T,D} <: StructuredGrid{T,D}  (AbstractRanges)
# RectilinearGrid{T,D} <: StructuredGrid{T,D} (sorted list for each Cartesian directions)
# struct CartesianGrid{Dim, T<:Real, V<:StaticVector{Dim,T}} <: AbstractGrid{Dim,V}
# end




using StaticArrays

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
struct CartesianGrid4{Dim,T<:AbstractFloat,V<:StaticVector{Dim,T}} <: AbstractArray{V,Dim}
    origin::V
    spacing::T
    dims::Dims{Dim}

    function CartesianGrid4{Dim,T,V}(origin, spacing, dims) where {Dim,T<:Real,V<:StaticVector{Dim,T}}
        @assert all(dims .> 0) "dimensions must be positive"
        @assert spacing > 0 "spacing must be positive"
        new(origin, spacing, dims)
    end
end

function CartesianGrid4(origin::StaticVector{N,T1}, spacing::T2, dims::Dims{N}) where {T1<:AbstractFloat,T2<:AbstractFloat,N}
    T = promote_type(T1, T2)
    CartesianGrid4(convert.(T, origin), convert(T, spacing), dims)
end

# AbstractVector but not StaticVector -> use Svector as fallback for default value
function CartesianGrid4(origin::AbstractVector{T1}, spacing::T2, dims::Dims{N}) where {T1<:AbstractFloat,T2<:AbstractFloat,N}
    @assert length(origin) == N "Length of origin must match number of dimensions"
    T = promote_type(T1, T2)
    CartesianGrid4(SVector{N,T}(origin...), convert(T, spacing), dims)
end

Base.size(g::CartesianGrid4) = g.dims
Base.ndims(::CartesianGrid4{Dim}) where {Dim} = Dim

Base.minimum(g::CartesianGrid4) = g.origin
Base.maximum(g::CartesianGrid4) = g.origin + width(g)
Base.extrema(g::CartesianGrid4) = (minimum(g), maximum(g))

spacing(g::CartesianGrid4) = g.spacing
width(g::CartesianGrid4{Dim,T,V}) where {Dim,T,V} = V((g.dims .- 1) * g.spacing)
center(g::CartesianGrid4) = g.origin + width(g) / 2

@inline function Base.getindex(g::CartesianGrid4{Dim}, I::Vararg{Int,Dim}) where {Dim}
    @boundscheck checkbounds(g, I...)
    return @. g.origin + (I - 1) * g.spacing
end

function Base.show(io::IO, g::CartesianGrid4{Dim,T,V}) where {Dim,T,V}
    dims = join(g.dims, "×")
    print(io, "$dims CartesianGrid4{$Dim,$T,$V}")
end

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid4)
    println(io, g)
    println(io, "  minimum: ", minimum(g))
    println(io, "  maximum: ", maximum(g))
    print(io, "  spacing: ", spacing(g))
end



g1 = CartesianGrid4(SVector{3,Float64}(1, 2, 3), 1, Dims((10, 10, 10)))
g2 = CartesianGrid4(SVector(1, 2, 3.0), 1.0, (10, 10, 10))
g3 = CartesianGrid4(SVector(1, 2, 3), 1.0, (10, 10, 10))  # missmatch type (promotion)
g4 = CartesianGrid4([1, 2, 3], 1.0, (10, 10, 10)) # Svector fallback as default
CartesianGrid4([1, 2, 3], 1.0, (3, 3, 3))




# # A regular grid is a tessellation of n-dimensional Euclidean space by congruent parallelotopes (e.g. bricks).
# # A Cartesian grid is a special case where the elements are unit squares or unit cubes, and the vertices are points on the integer lattice.

# TODO: origin(start) spacing (steps) and dims (length)  VS start stop and length (dims)
# abstract type AbstractGrid{Dim,V} end
# abstract type StructuredGrid{Dim,V} <: AbstractGrid{Dim,V} end
# abstract type UnstructuredGrid{Dim,V} <: AbstractGrid{Dim,V} end

# CartesianGrid{T,D} <: StructuredGrid{T,D} (spacing, dims)
# RegularGrid{T,D} <: StructuredGrid{T,D}  (AbstractRanges)
# RectilinearGrid{T,D} <: StructuredGrid{T,D} (sorted list for each Cartesian directions)
# struct CartesianGrid{Dim, T<:Real, V<:StaticVector{Dim,T}} <: AbstractGrid{Dim,V}
# end





using StaticArrays

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
struct CartesianGrid2{Dim,T<:Real,V<:StaticVector{Dim,T}} <: AbstractArray{V,Dim}
    origin::V
    spacing::T
    dims::Dims{Dim}

    function CartesianGrid2{Dim,T,V}(origin, spacing, dims) where {Dim,T<:Real,V<:StaticVector{Dim,T}}
        @assert all(dims .> 0) "dimensions must be positive"
        @assert spacing > 0 "spacing must be positive"
        new(origin, spacing, dims)
    end
end

function CartesianGrid2(origin::StaticVector{N,T1}, spacing::T2, dims::Dims{N}) where {T1<:Real,T2<:Real,N}
    T = promote_type(T1, T2)
    CartesianGrid2(convert.(T, origin), convert(T, spacing), dims)
end

# AbstractVector but not StaticVector -> use Svector as fallback for default value
function CartesianGrid2(origin::AbstractVector{T1}, spacing::T2, dims::Dims{N}) where {T1<:Real,T2<:Real,N}
    @assert length(origin) == N "Length of origin must match number of dimensions"
    T = promote_type(T1, T2)
    CartesianGrid2(SVector{N,T}(origin...), convert(T, spacing), dims)
end

Base.size(g::CartesianGrid2) = g.dims
Base.ndims(::CartesianGrid2{Dim}) where {Dim} = Dim

Base.minimum(g::CartesianGrid2) = g.origin
Base.maximum(g::CartesianGrid2) = g.origin + width(g)
Base.extrema(g::CartesianGrid2) = (minimum(g), maximum(g))

spacing(g::CartesianGrid2) = g.spacing
width(g::CartesianGrid2{Dim,T,V}) where {Dim,T,V} = V((g.dims .- 1) .* g.spacing)
center(g::CartesianGrid2) = g.origin + width(g) ./ 2

@inline function Base.getindex(g::CartesianGrid2{Dim}, I::Vararg{Int,Dim}) where {Dim}
    @boundscheck checkbounds(g, I...)
    return @. g.origin + (I - 1) * g.spacing
end

function Base.show(io::IO, g::CartesianGrid2{Dim,T,V}) where {Dim,T,V}
    dims = join(g.dims, "×")
    print(io, "$dims CartesianGrid2{$Dim,$T,$V}")
end

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid2)
    println(io, g)
    println(io, "  minimum: ", minimum(g))
    println(io, "  maximum: ", maximum(g))
    print(io, "  spacing: ", spacing(g))
end



g1 = CartesianGrid2(SVector{3,Int64}(1, 2, 3), 1, Dims((10, 10, 10)))
g2 = CartesianGrid2(SVector(1, 2, 3.0), 1.0, (10, 10, 10))
g3 = CartesianGrid2(SVector(1, 2, 3), 1.0, (10, 10, 10))  # missmatch type (promotion)
g4 = CartesianGrid2([1, 2, 3], 1.0, (10, 10, 10)) # Svector fallback as default
CartesianGrid2([1, 2, 3], 1.0, (3, 3, 3))




# # A regular grid is a tessellation of n-dimensional Euclidean space by congruent parallelotopes (e.g. bricks).
# # A Cartesian grid is a special case where the elements are unit squares or unit cubes, and the vertices are points on the integer lattice.

# TODO: origin(start) spacing (steps) and dims (length)  VS start stop and length (dims)
# abstract type AbstractGrid{Dim,V} end
# abstract type StructuredGrid{Dim,V} <: AbstractGrid{Dim,V} end
# abstract type UnstructuredGrid{Dim,V} <: AbstractGrid{Dim,V} end

# CartesianGrid{T,D} <: StructuredGrid{T,D} (spacing, dims)
# RegularGrid{T,D} <: StructuredGrid{T,D}  (AbstractRanges)
# RectilinearGrid{T,D} <: StructuredGrid{T,D} (sorted list for each Cartesian directions)
# struct CartesianGrid{Dim, T<:Real, V<:StaticVector{Dim,T}} <: AbstractGrid{Dim,V}
# end
