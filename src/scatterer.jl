abstract type Scatter end

# TODO:
# - defintion of material of the scatterers

struct Sphere <: Scatter
    radius::Float64
    center::SVector{3,Float64}
end


function isinside(x::Real, y::Real, z::Real, s::Sphere)
    r = s.radius
    x0, y0, z0 = s.center
    if abs2(x - x0) + abs2(y - y0) + abs2(z - z0) <= abs2(r)
        return true
    end
    return false
end


struct Disk <: Scatter
    radius::Float64
    height::Float64
    center::SVector{3,Float64}
end


function isinside(x::Real, y::Real, z::Real, s::Disk)
    r, h = s.radius, s.height
    x0, y0, z0 = s.center
    if abs2(x - x0) + abs2(y - y0) <= abs2(r) &&
        z >= z0 && z <= z0 + h
        return true
    end
    return false
end


function dipoles(g::CartesianGrid, s::Scatter)
    inds = CartesianIndices(g)
    out = SVector{3,Float64}[]
    for ind in inds
        coord = g[ind]
        if isinside(coord..., s)
            push!(out, coord)
        end
    end
    return out
end



function discretize(g::CartesianGrid, s::Scatter)

    eps = 1
    out = zeros(typeof(eps), size(g))

    inds = CartesianIndices(g)
    for ind in inds
        coord = g[ind]
        if isinside(coord..., s)
            out[ind] = eps
        end
    end
    return out
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












