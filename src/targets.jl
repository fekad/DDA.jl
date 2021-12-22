abstract type AbstractShape end

struct Dipole <: AbstractShape
    center::Point3{Float64}
end

function dipoles(g::CartesianGrid{T,N}, s::Dipole) where {T,N}

    ind = @. round(Int, (s.center - g.origin) / g.spacing + 1)
    inds = [CartesianIndex(Tuple(ind)), ]
    coords = g[inds]

    return coords, inds
end


struct Sphere <: AbstractShape
    center::Point3{Float64}
    radius::Float64
end

function Base.in(p::Point, s::Sphere)
    return norm(s.center - p) ≤ s.radius
end


struct Disk <: AbstractShape
    center::Point3{Float64}
    radius::Float64
    height::Float64
end


function Base.in(p::Point, s::Disk)
    x, y, z = p
    r, h = s.radius, s.height
    x0, y0, z0 = s.center
    return  abs2(x - x0) + abs2(y - y0) <= abs2(r) && z >= z0 && z <= z0 + h
end


# TODO: adding center
struct TriangularPlatelets <: AbstractShape
    width::Float64
    height::Float64
end

function Base.in(p::Point, s::TriangularPlatelets)
    x, y, z = p
    w, h = s.width, s.height
    return  x >= - sqrt(3) / 4 * w && abs(y) <= w/4 - x/sqrt(3) && z >= 0 && z <= h
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
    return coords,inds
end



function discretize(g::CartesianGrid, s::AbstractShape)

    occupation = zeros(Bool, size(g))

    for i in eachindex(g)
        coordinate = g[i]

        if coordinate ∈ s
            occupation[i] = true
        end
    end
    return occupation
end


# function dipoles(g::CubicGrid, s::Sphere)
#     return [[x, y, z] for x in g.xrange, y in g.yrange, z in g.zrange if sqrt(x^2 + y^2 + z^2) <= s.radius]
# end
#
#
# function discretize(g::CubicGrid, s::Sphere, eps)
#     out = zeros(typeof(eps), dims(g))
#     for (i, x) in enumerate(g.xrange),
#         (j, y) in enumerate(g.yrange),
#         (k, z) in enumerate(g.zrange)
#         # @show x, xi
#         if sqrt(x^2 + y^2 + z^2) <= s.radius
#             out[i,j,k] = eps
#         end
#     end
#     return out
# end
#
#

# g = CubicGrid(-1, 1, -1, 1, -1, 1, .1)
# s = Sphere([0,0,0], .5)

# dipoles = discretize(g, s, 1.0 + .0im)
# @time pos = positions(g, s)










# function isinside(p::Point, s::Sphere)
#     # x0, y0, z0 = s.center
#     # x, y, z = p
#     # return abs2(x - x0) + abs2(y - y0) + abs2(z - z0) <= abs2(s.radius)
#     # return sum(@. abs2(s.center - p)) ≤ abs2(s.radius)
#     return norm(s.center - p) ≤ s.radius
# end

# function isinside(x::Real, y::Real, z::Real, s::Sphere)
#     r = s.radius
#     x0, y0, z0 = s.center
#     return abs2(x - x0) + abs2(y - y0) + abs2(z - z0) <= abs2(r)
# end



# function dipoles(g::CartesianGrid{T,N}, s::AbstractShape) where {T,N}
#     coords = Point{N,T}[]
#     for ind in CartesianIndices(g)
#         coord = g[ind]
#         if isinside(coord..., s)
#             push!(coords, coord)
#         end
#     end
#     return coords
# end

# function dipoles2(g::CartesianGrid{T,N}, s::AbstractShape) where {T,N}
#
#     coords = Point{N,T}[]
#     inds = CartesianIndex{N}[]
#
#     for ind in CartesianIndices(g)
#         coord = g[ind]
#         if isinside(coord..., s)
#             push!(coords, coord)
#             push!(inds, ind)
#         end
#     end
#     return coords,inds
# end
















# sphere(x::Real, y::Real, z::Real, r::Real) = abs2(x) + abs2(y) + abs2(z) <= abs2(r)
# sphere(x::Real, y::Real, z::Real, s::Sphere) = sphere(x - s.center[1], y - s.center[2], z - s.center[3], s.radius)
#
# function dipoles(g::CartesianGrid, s::Sphere)
#     inds = CartesianIndices(g)
#     return [g[ind] for ind in inds if sphere(g[ind]..., s)]
# end





# function dipoles(g::CubicGrid, s::Sphere, eps)
#     inds = CartesianIndices(g)
#     out = zeros(typeof(eps), dims(g))
#     @inbounds for i in inds
#         # @show x, xi
#         x, y, z = g[i]
#         if sqrt(x^2 + y^2 + z^2) <= s.radius
#             out[i,j,k] = eps
#         end
#     end
#     return out
# end

# g = CartesianGrid([.5,.5,.5], ones(3), (30,30,30))
# s = Sphere(15., [15,15,15])
# d = dipoles(g,s)
#
#
# g = CartesianGrid([-14.5, -14.5, -14.5], ones(3), (30,30,30))
# s = Sphere(15., [0, 0, 0])
# d = dipoles(g,s)
#
#
# scatter3d(dd[1,:],dd[2,:],dd[3,:])
#












