# Working but it could b e nicer
#using AbstractPlotting
using Colors
using GeometryBasics

using WGLMakie
WGLMakie.activate!()

# using GLMakie
# GLMakie.inline!(false)
# GLMakie.activate!()
# f = Figure();
# sc = display(f);
# ax = Axis3(f[1, 1], title = "Figure 1", aspect=(1, 1, 1), perspectiveness=0.);
# hidedecorations!(ax)

scene = Scene(camera = cam3d!)

# Figure options
# fig_size = [500, 500]
# fig_view = [45, 65] # Azimuthal and Elvation viewing angles

# Sphere options
sphere_color = (:grey80, .6)

arrowsize=0.2

# Wireframe options
frame_num = 4
frame_color =  (:black, .2)
frame_radius = 1.

# Axes
axes_color = (:black, .4)
axes_radius = 3.0


sphere = Tesselation(Sphere(Point3f0(0), 1f0), 100)
sphere_mesh = uv_normal_mesh(sphere)


mesh!(scene, sphere_mesh, color=sphere_color, show_axis = false)


# θ ∈ [0, π],  ϕ ∈ [0, 2π]
inner(θ, ϕ) = cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)

function latcurves(θ)
    nvertices = 100
    ϕ = LinRange(0, 2π, nvertices)
    return [inner(θ, ϕ) for θ in θ for θ in θ, ϕ in ϕ ]
end

function loncurves(ϕ)
    nvertices = 100
    θ = LinRange(0, π, nvertices)
    return [inner(θ, ϕ) for θ in θ for θ in θ, ϕ in ϕ ]
end


for θ = LinRange(0, π, 2*frame_num+1)
    lines!(scene, latcurves(θ),
        color = frame_color,
        linewidth = frame_radius
    )
end
for ϕ = LinRange(0, 2π, 2*frame_num+1)
    lines!(scene, loncurves(ϕ),
        color = frame_color,
        linewidth = frame_radius
    )
end


# add axes
axis = [-1, 1]
other = [0, 0]
lines!(scene,
    axis, other, other,
    color = axes_color,
    linewidth = axes_radius
)

lines!(scene,
    other, axis, other,
    color = axes_color,
    linewidth = axes_radius
)
lines!(scene,
    other, other, axis,
    color = axes_color,
    linewidth = axes_radius
)


# Labels
text!(scene, "|R>", textsize=30, position = Point3f0(0, 0,  1.4), align = (:center, :center))
text!(scene, "|L>", textsize=30, position = Point3f0(0, 0, -1.4), align = (:center, :center))

text!(scene, "|H>", textsize=30, position = Point3f0( 1.4, 0, 0), align = (:center, :center))
text!(scene, "|V>", textsize=30, position = Point3f0(-1.4, 0, 0), align = (:center, :center))

text!(scene, "|D>", textsize=30, position = Point3f0(0, 1.4, 0), align = (:center, :center))
text!(scene, "|A>", textsize=30, position = Point3f0(0,-1.4, 0), align = (:center, :center))



# arrows!(scene, [Point3f0(0)], [Point3f0(0,0,1 - arrowsize)], arrowsize=arrowsize)
arrows!(scene, [Point3f0(0)], [Point3f0(1 - arrowsize,0,0)], arrowsize=arrowsize, color = :red)

display(scene)

