# using Base: Float64, disable_library_threading_hooks
using LazyGrids: ndgrid

abstract type AbstractGrid end

import Base:length

struct CubicGrid <: AbstractGrid
    xrange::AbstractRange
    yrange::AbstractRange
    zrange::AbstractRange
    function CubicGrid(xrange, yrange, zrange)
        @assert isapprox(step(xrange), step(yrange))
        @assert isapprox(step(xrange), step(zrange))
        return new(xrange, yrange, zrange)
    end
end

CubicGrid(xmin, xmax, ymin, ymax, zmin, zmax, dx) = CubicGrid(range(xmin, xmax, step=dx), range(ymin, ymax, step=dx), range(zmin, zmax, step=dx))

dims(g::CubicGrid) = (length(g.xrange), length(g.yrange), length(g.zrange))
length(g::CubicGrid) = *(dims(g)...)


abstract type Scatter end

struct Sphere <: Scatter
    radius::Float64
    center::Vector
end

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
eltype(CubicGrid)

function positions(g::CubicGrid, s::Sphere)
    out = Array{Float64}(undef, (3, length(g)))

    # for i eachindex(g)

    return [[x, y, z] for x in g.xrange, y in g.yrange, z in g.zrange if sqrt(x^2 + y^2 + z^2) <= s.radius]

end

g = CubicGrid(-1, 1, -1, 1, -1, 1, .1)
# s = Sphere([0,0,0], .5)

# dipoles = discretize(g, s, 1.0 + .0im)
# @time pos = positions(g, s)



struct Disk <: Scatter
    radius::Float64
    height::Float64
    center::Vector
    orientation::Vector
end







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
#
#
# struct CubicGrid
#     step
#     xmin
#     xmax
#     ymin
#     ymax
#     zmin
# end
#
#
#
#
# function get_coords(grid::Grid)
#
# end


#

#
# x = Linrange(0, 10, 1001)
#
# CubicLattice(x, y, z)