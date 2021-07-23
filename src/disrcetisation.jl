# using Base: Float64, disable_library_threading_hooks

using StaticArrays

abstract type AbstractGrid end

struct CubicGrid <: AbstractGrid
    xi::StepRangeLen
    yi::StepRangeLen
    zi::StepRangeLen
end

origin(g::CubicGrid) = [g.xi[1], g.yi[1], g.zi[1]]
dims(g::CubicGrid) = (length(g.xi), length(g.yi), length(g.zi))

abstract type Scatter end

struct Sphere <: Scatter
    grid::CubicGrid
    center::Vector
    radius::Float64
end

function grid(s::Sphere)
    out = zeros(Float64, dims(s.grid))
    for (i, x) in enumerate(s.grid.xi), (j, y) in enumerate(s.grid.yi), (k, z) in enumerate(s.grid.zi)
        # @show x, xi
        if sqrt(x^2 + y^2 + z^2) <= s.radius
            out[i,j,k] = 1.
        end
    end
    return out
end

function dipoles(s::Sphere)

end

# cubic grid
# - step size (dx == dy == dz)

# - start(3) stop(3) step()
# - start(3) stop(3)

struct CubicLattice <: AbstractGrid
    origin::SVector{3,Float64}
    dims::SVector{3,Int64}
    spacing::Float64
end

struct RotatadCubicLattice <: AbstractGrid
    start::SVector{3,Float64}
    length::SVector{3,Int64}
    step::Float64
    orientation::SMatrix{3,3,Float64}
end



struct CubicGrid
    step
    xmin
    xmax
    ymin
    ymax
    zmin
end




function get_coords(grid::Grid)

end


sphere
disk



x = Linrange(0, 10, 1001)

CubicLattice(x, y, z)