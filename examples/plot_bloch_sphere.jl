using Plots
plotlyjs()
# gr()

# using GeometryBasics
# s = GeometryBasics.Sphere(Point3(0),1)
# p = GeometryBasics.coordinates(s, 24)


n = 21
ϕ = range(-π, π; length = n)
θ = range(0, π; length = n)
x = cos.(ϕ) * sin.(θ)'
y = sin.(ϕ) * sin.(θ)'
z = ones(n) * cos.(θ)'

# wireframe(x, y, z)
plot(x, y, z; seriestype=:surface)
plot(x, y, z; seriestype=:wireframe)

