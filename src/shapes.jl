abstract type AbstractShape{Dim,T} end

struct Dipole{Dim,T} <: AbstractShape{Dim,T}
    origin::StaticVector{Dim,T}
end

function dipoles(g::CartesianGrid{Dim,T}, s::Dipole{Dim,T}) where {Dim,T}

    ind = @. round(Int, (s.origin - g.origin) / g.spacing + 1)
    inds = [CartesianIndex(Tuple(ind)),]
    coords = g[inds]

    return coords, inds
end


struct Sphere{Dim,T} <: AbstractShape{Dim,T}
    origin::StaticVector{Dim,T}
    radius::T
end

function Base.in(p::StaticVector{Dim,T}, s::Sphere{Dim,T}) where {Dim,T}
    return norm(p - s.origin) ≤ s.radius
end


struct Disk{Dim,T} <: AbstractShape{Dim,T}
    origin::StaticVector{Dim,T}
    radius::T
    height::T
end


function Base.in(p::StaticVector{Dim,T}, s::Disk{Dim,T}) where {Dim,T}
    x, y, z = p - s.origin
    r, h = s.radius, s.height
    return abs2(x) + abs2(y) <= abs2(r) && z >= -h / 2 && z <= h / 2
end


struct Box{Dim,T} <: AbstractShape{Dim,T}
    origin::StaticVector{Dim,T}
    size::SVector{Dim,T}
end

function Base.in(p::StaticVector{Dim,T}, s::Box{Dim,T}) where {Dim,T}
    return all(@. 0 ≤ p - s.origin ≤ s.size)
end


struct TriangularPlatelets{Dim,T} <: AbstractShape{Dim,T}
    origin::StaticVector{Dim,T}
    width::T
    height::T
end

function Base.in(p::StaticVector{Dim,T}, s::TriangularPlatelets{Dim,T}) where {Dim,T}
    x, y, z = p - s.origin
    w, h = s.width, s.height
    return x >= -sqrt(3) / 4 * w && abs(y) <= w / 4 - x / sqrt(3) && z >= 0 && z <= h
end


struct Composite{Dim,T,I<:AbstractShape{Dim,T}} <: AbstractShape{Dim,T}
    shapes::Vector{I}
end

Base.getindex(multi::Composite, ind) = getindex(multi.shapes, ind)
Base.length(multi::Composite) = length(multi.shapes)
Base.eltype(multi::Composite) = eltype(multi.shapes)
Base.firstindex(multi::Composite) = firstindex(multi.shapes)
Base.lastindex(multi::Composite) = lastindex(multi.shapes)
Base.iterate(multi::Composite, state = 1) =
    state > length(multi) ? nothing : (multi[state], state + 1)


function Base.in(p::StaticVector{Dim,T}, s::Composite{Dim,T}) where {Dim,T}
    return mapreduce(s -> in(p, s), |, s)
end


function dipoles(g::CartesianGrid{Dim,T,V}, s::AbstractShape{Dim,T}) where {Dim,T,V}

    coords = V[]

    for i in eachindex(g)
        coordinate = g[i]
        if coordinate ∈ s
            push!(coords, coordinate)
        end
    end
    return coords
end

function indices(g::CartesianGrid{Dim,T,V}, s::AbstractShape{Dim,T}) where {Dim,T,V}

    inds = CartesianIndex{Dim}[]

    for i in eachindex(g)
        coordinate = g[i]
        if coordinate ∈ s
            push!(inds, i)
        end
    end
    return inds
end

function discretize(g::CartesianGrid, s::S) where {S<:AbstractShape}

    occupation = zeros(Bool, size(g))

    for i in eachindex(g)
        coordinate = g[i]
        if coordinate ∈ s
            occupation[i] = true
        end
    end
    return occupation
end

