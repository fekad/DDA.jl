

abstract type Scatter end

struct Sphere <: Scatter
    radius::Float64
    center::Vector
end

sphere(x::Real, y::Real, z::Real, r::Real) = abs2(x) + abs2(y) + abs2(z) <= abs2(r)
sphere(x::Real, y::Real, z::Real, s::Sphere) = sphere(x - s.center[1], y - s.center[2], z - s.center[3], s.radius)

function dipoles(g::CartesianGrid, s::Sphere)
    inds = CartesianIndices(g)
    return [g[ind] for ind in inds if sphere(g[ind]..., s)]
end


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

















# sphere
# disk

function discretize(g::CubicGrid, s::Sphere, eps)
    out = zeros(typeof(eps), dims(g))
    for (i, x) in enumerate(g.xrange),
        (j, y) in enumerate(g.yrange),
        (k, z) in enumerate(g.zrange)
        # @show x, xi
        if sqrt(x^2 + y^2 + z^2) <= s.radius
            out[i,j,k] = eps
        end
    end
    return out
end
# eltype(CubicGrid)

function positions(g::CubicGrid, s::Sphere)
    out = Array{Float64}(undef, (3, length(g)))

    # for i eachindex(g)

    return [[x, y, z] for x in g.xrange, y in g.yrange, z in g.zrange if sqrt(x^2 + y^2 + z^2) <= s.radius]

end

# g = CubicGrid(-1, 1, -1, 1, -1, 1, .1)
# s = Sphere([0,0,0], .5)

# dipoles = discretize(g, s, 1.0 + .0im)
# @time pos = positions(g, s)



struct Disk <: Scatter
    radius::Float64
    height::Float64
    center::Vector
    orientation::Vector
end
