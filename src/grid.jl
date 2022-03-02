
"""
    CartesianGrid(origin, spacing, dims)

A Cartesian grid with dimensions `dims`, lower left corner at `origin`
and cell spacing `spacing`.

## Examples

Create a 2D grid with 100x100 locations and origin at [10.,20.] units:

```julia
julia> CartesianGrid([10.,20.], 1., (100, 100))
```
"""
struct CartesianGrid{Dim,T<:Real,V<:StaticVector{Dim,T}} <: AbstractArray{V,Dim}
    origin::V
    spacing::T
    dims::Dims{Dim}

    function CartesianGrid{Dim,T,V}(origin, spacing, dims) where {Dim,T<:Real,V<:StaticVector{Dim,T}}
        @assert all(dims .> 0) "dimensions must be positive"
        @assert spacing > 0 "spacing must be positive"
        new(origin, spacing, dims)
    end
end

CartesianGrid(origin::V, spacing::T, dims::Dims{Dim}) where {Dim,T<:Real,V<:StaticVector{Dim,T}} = CartesianGrid{Dim,T,V}(origin, spacing, dims)

function CartesianGrid(origin::StaticVector{N,T1}, spacing::T2, dims::Dims{N}) where {T1<:Real,T2<:Real,N}
    T = promote_type(T1, T2)
    CartesianGrid(convert.(T, origin), convert(T, spacing), dims)
end

# AbstractVector but not StaticVector -> use Svector as fallback for default value
function CartesianGrid(origin::AbstractVector{T1}, spacing::T2, dims::Dims{N}) where {T1<:Real,T2<:Real,N}
    @assert length(origin) == N "Length of origin must match number of dimensions"
    T = promote_type(T1, T2)
    CartesianGrid(SVector{N,T}(origin...), convert(T, spacing), dims)
end

Base.size(g::CartesianGrid) = g.dims
Base.ndims(::CartesianGrid{Dim}) where {Dim} = Dim

Base.minimum(g::CartesianGrid) = g.origin
Base.maximum(g::CartesianGrid) = g.origin + width(g)
Base.extrema(g::CartesianGrid) = (minimum(g), maximum(g))

spacing(g::CartesianGrid) = g.spacing
width(g::CartesianGrid{Dim,T,V}) where {Dim,T,V} = V((g.dims .- 1) .* g.spacing)
center(g::CartesianGrid) = g.origin + width(g) ./ 2

@inline function Base.getindex(g::CartesianGrid{Dim}, I::Vararg{Int,Dim}) where {Dim}
    @boundscheck checkbounds(g, I...)
    return @. g.origin + (I - 1) * g.spacing
end

function Base.show(io::IO, g::CartesianGrid{Dim,T,V}) where {Dim,T,V}
    dims = join(g.dims, "×")
    print(io, "$dims CartesianGrid{$Dim,$T,$V}")
end

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid)
    println(io, g)
    println(io, "  minimum: ", minimum(g))
    println(io, "  maximum: ", maximum(g))
    print(io, "  spacing: ", spacing(g))
end


# CartesianGrid(SVector{3,Int64}(1, 2, 3), 1, Dims((10, 10, 10)))
# CartesianGrid(SVector(1, 2, 3.0), 1.0, (10, 10, 10))
# CartesianGrid(SVector(1, 2, 3), 1.0, (10, 10, 10))  # missmatch type (promotion)
# CartesianGrid([1, 2, 3], 1.0, (10, 10, 10)) # Svector fallback as default
