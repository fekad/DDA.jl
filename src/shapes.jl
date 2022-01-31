
struct Dipole <: AbstractShape
    origin::Point3{Float64}
end

function dipoles(g::CartesianGrid{T,N}, s::Dipole) where {T,N}

    ind = @. round(Int, (s.origin - g.origin) / g.spacing + 1)
    inds = [CartesianIndex(Tuple(ind)),]
    coords = g[inds]

    return coords, inds
end


struct Sphere <: AbstractShape
    origin::Point3{Float64}
    radius::Float64
end

function Base.in(p::Point, s::Sphere)
    return norm(p - s.origin) ≤ s.radius
end


struct Disk <: AbstractShape
    origin::Point3{Float64}
    radius::Float64
    height::Float64
end


function Base.in(p::Point, s::Disk)
    x, y, z = p - s.origin
    r, h = s.radius, s.height
    return abs2(x) + abs2(y) <= abs2(r) && z >= -h/2 && z <= h/2
end


struct Box <: AbstractShape
    origin::Point3{Float64}
    size::SVector{3,Float64}
end

function Base.in(p::Point, s::Box)
    return all(@. 0 ≤ p - s.origin ≤ s.size)
end


struct TriangularPlatelets <: AbstractShape
    origin::Point3{Float64}
    width::Float64
    height::Float64
end

function Base.in(p::Point, s::TriangularPlatelets)
    x, y, z = p - s.origin
    w, h = s.width, s.height
    return x >= -sqrt(3) / 4 * w && abs(y) <= w / 4 - x / sqrt(3) && z >= 0 && z <= h
end


struct Composite{I<:AbstractShape} <:AbstractShape
    shapes::Vector{I}
end

Base.getindex(multi::Composite, ind) = getindex(multi.shapes, ind)
Base.length(multi::Composite) = length(multi.shapes)
Base.eltype(multi::Composite) = eltype(multi.shapes)
Base.firstindex(multi::Composite) = firstindex(multi.shapes)
Base.lastindex(multi::Composite) = lastindex(multi.shapes)
Base.iterate(multi::Composite, state=1) =
  state > length(multi) ? nothing : (multi[state], state+1)


function Base.in(p::Point, s::Composite)
    return mapreduce(s -> in(p, s), |, s)
end



function dipoles(g::CartesianGrid{T,N}, s::AbstractShape) where {T,N}

    coords = Point{N,T}[]
    inds = CartesianIndex{N}[]

    for ind in CartesianIndices(g)
        coord = g[ind]
        if coord ∈ s
            push!(coords, coord)
            push!(inds, ind)
        end
    end
    return coords, inds
end



function discretize(g::CartesianGrid, s::S) where S<:AbstractShape

    occupation = zeros(Bool, size(g))

    for i in eachindex(g)
        coordinate = g[i]

        if coordinate ∈ s
            occupation[i] = true
        end
    end
    return occupation
end

