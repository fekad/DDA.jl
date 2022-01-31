

using DelimitedFiles: readdlm

filename = "examples/triangle.pos"

coords = readdlm(filename, comments=true)[:,1:3]

using  Plots
plotlyjs()


scatter(eachcol(coords)...)

using DDA
using StaticArrays

# 1686
# Triangular platelets (lateral size = 10 nm, thickness = 2.4 nm). 6 layers of 281 dipoles
# This should be 7 layers

origin = Point3(-48., -48, 0)
spacing = [4., 4., 4.]
dims=(25, 25, 7)
grid = DDA.CartesianGrid(origin, spacing, dims)




# Equilateral Triangle
# struct TriangularPlatelets <: DDA.AbstractTarget
#     center::Point3{Float64}
#     width::Float64
#     height::Float64
# end
struct TriangularPlatelets2 <: DDA.AbstractTarget
    width::Float64
    height::Float64
end

A_ddeels = 1686/6 *4^2
A_analytic = sqrt(3)/4 * 100^2

function Base.in(p::Point, s::TriangularPlatelets2)
    x, y, z = p
    w, h = s.width, s.height
    # return  x >= - sqrt(3) / 4 * w && abs(y) - .6 <= w/4 - x/sqrt(3) && z >= 0 && z <= h
    return  x >= - sqrt(3) / 4 * w && abs(y) <= w/4 - x/sqrt(3) && z >= 0 && z <= h
end

s = TriangularPlatelets2(100., 1.)
dipoles, _ = DDA.dipoles(grid, s)


Nd = length(dipoles) #263
A = Nd * 4^2

A_ddeels/A_analytic*100, A/A_analytic*100

x = [d[1] for d in dipoles]
y = [d[2] for d in dipoles]
scatter(x, y ; markersize=6, aspect_ratio=:equal, label="generated");

# xx,yy = x,y
# scatter!(xx,yy; markersize = 4, label="ddeels example")
scatter!(coords[:,1], coords[:,2]; markersize = 4, label="ddeels example")

xr = -48:.1:48
f(x) = 100/4 - x/sqrt(3) + 2
fig=plot!(xr,f.(xr))


savefig(fig, "examples/triangle_plus_2.pdf")








