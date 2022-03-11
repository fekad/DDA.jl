
using DDA

using FFTW
# FFTW.set_provider!("fftw")
FFTW.set_num_threads(1)

using StaticArrays
# using GeometryBasics

using Unitful
using PhysicalConstants.CODATA2018: h, c_0
using RefractiveIndexDatabase
using Permittivity

using StaticArrays
using IterativeSolvers
using LinearAlgebra: norm

using TimerOutputs

using PlotlyJS
# using Plots
# plotlyjs()
# # gr()

# Au = get_material("main", "Au", "Johnson")
Au = get_material("main", "Au", "Babar")
# Au = get_material("main", "Au", "Olmon-sc")
# Au = get_material("main", "Au", "Olmon-ev")

# Convert into permititvity

λ_exp = Au.λ * 1u"μm" .|> u"nm"
m_exp = Au.n + Au.k * im
ε_exp = m_exp .^ 2

# interpolation
model = PermittivityTable(ustrip.(λ_exp), ε_exp)

##

l = range(620, 740, step = 1.0)

# Drude model
# DOI: 10.1002/adom.201300457
ε_inf = 2.2715
ω_p = 8.9234u"eV"
Γ_p = 0.042389u"eV"

# optimised
ε_inf = 9.751684946459113
ω_p = 8.955375828538855u"eV"
Γ_p = 0.05151667730721536u"eV"

model_drude = DrudeModel(ε_inf, ustrip(ω_p), ustrip(Γ_p))
ω = h * c_0 ./ (l * u"nm") .|> u"eV"
ε_drude = model_drude.(ustrip(ω))




trace = [
    scatter(x = l, y = real(model(l)), mode="markers"),
    scatter(x = l, y = imag(model(l)), mode="markers", yaxis = "y2"),
    scatter(x = l, y = real(ε_drude), mode="markers"),
    scatter(x = l, y = imag(ε_drude), mode="markers", yaxis = "y2")
]

layout = Layout(
    title = "Delectric Permittivity",
    xaxis = attr(title = "Wavelength (nm)"),
    yaxis = attr(
        title = "real",
        titlefont = attr(color = :blue),
        tickfont = attr(color = :blue),
    ),
    yaxis2 = attr(
        title = "imag",
        titlefont = attr(color = :red),
        tickfont = attr(color = :red),
        overlaying = "y",
        side = :right
    )
)
plot(trace, layout)



## C_ext


# 1. Define a grid
disk_d = 10.0 # nm
disk_h = 1.0 # nm
disk_gap = 4 / 3 # nm

spacing = 1 / 3

Nx = round(Int, (2 * disk_d + disk_gap) / spacing)
Ny = round(Int, disk_d / spacing)
Nz = round(Int, disk_h / spacing)

# offset = [spacing / 2, spacing / 2, spacing / 2]
# origin = [disk_d / 2, disk_d / 2 , disk_h / 2]
origin = [0, 0, 0]

grid = CartesianGrid(origin, spacing, (Nx, Ny, Nz))
grid


# 2. Define the target(s)
# offset = [spacing / 2, spacing / 2, spacing / 2]
radius = disk_d / 2

origin = SVector(disk_d / 2 - spacing / 2, disk_d / 2 - spacing / 2, disk_h / 2)
disk1 = DDA.Disk(origin, disk_d / 2, disk_h)

origin = SVector(disk_d + disk_gap + disk_d / 2 - spacing / 2, disk_d / 2 - spacing / 2, disk_h / 2)
disk2 = DDA.Disk(origin, disk_d / 2, disk_h)

inds1 = indices(grid, disk1)
inds2 = indices(grid, disk2)


inds = [inds1..., inds2...]
coords = grid[inds]

layout = Layout(scene = attr(aspectmode = :data, camera_projection_type = :orthographic))
trace = scatter3d(x = [c[1] for c in coords], y = [c[2] for c in coords], z = [c[3] for c in coords], marker = attr(color = 1:length(coords)), mode = "markers")
plot(trace, layout)

# Parameters
e = [1, 0]        # Jones polarisation vector
θ, ϕ = 0.0, 0.0   # rotation angles [rad]

lambda = 600 # nm
k = 2π / lambda
ε = model(lambda)

kvec = [0, 0, k]
E₀ = [e..., 0]
alpha = LDR(ε, spacing, kvec, E₀)
alphas = fill(alpha, size(inds))


# 3. calculated the incident field Einc, at each dipole,
pw = PlaneWave(k, e, θ, ϕ)
Eincs = [field(pw, coord) for coord in coords]


layout = Layout(scene = attr(aspectmode = :data, camera_projection_type = :orthographic))
trace = cone(
    x = [c[1] for c in coords],
    y = [c[2] for c in coords],
    z = [c[3] for c in coords],
    u = [real(c[1]) for c in Eincs],
    v = [real(c[2]) for c in Eincs],
    w = [real(c[3]) for c in Eincs],
)
plot(trace, layout)


# 5. Define the DDA problem
prob = GridProblem(k, norm(e), grid, inds, alphas, Eincs)

# 6. Solve the DDA problem
sol = solve(prob, BiCGStabl(), reltol = 1e-8, verbose = true)

@show  C_abs(sol)
@show  C_sca(sol)


##

# quiver(Tuple.(coords), quiver = Tuple.(real.(Eincs)), aspect_ratio = :equal)

# Tuple.(real.(Eincs))


# Parameters:
lambda_range = 632:1.0:732 # nm

Q_abs = zeros(length(lambda_range))
Q_sca = zeros(length(lambda_range))


function build_system()


    # 1. Define a grid
    disk_d = 10.0 # nm
    disk_h = 1.0 # nm
    disk_gap = 20 / 3 # nm

    spacing = 1 / 3

    Nx = round(Int, (2 * disk_d + disk_gap) / spacing)
    Ny = round(Int, disk_d / spacing)
    Nz = round(Int, disk_h / spacing)

    # offset = [spacing / 2, spacing / 2, spacing / 2]
    # origin = [disk_d / 2, disk_d / 2 , disk_h / 2]
    origin = [0, 0, 0]

    grid = CartesianGrid(origin, spacing, (Nx, Ny, Nz))
    grid


    # 2. Define the target(s)
    # offset = [spacing / 2, spacing / 2, spacing / 2]
    radius = disk_d / 2

    origin = SVector(disk_d / 2 - spacing / 2, disk_d / 2 - spacing / 2, disk_h / 2)
    disk1 = DDA.Disk(origin, disk_d / 2, disk_h)

    origin = SVector(disk_d + disk_gap + disk_d / 2 - spacing / 2, disk_d / 2 - spacing / 2, disk_h / 2)
    disk2 = DDA.Disk(origin, disk_d / 2, disk_h)

    inds1 = indices(grid, disk1)
    inds2 = indices(grid, disk2)

    inds = [inds1..., inds2...]
    coords = grid[inds]

    return grid, inds, coords
end


function run_calc(lambda, model, verbose = false)
    grid, inds, coords = build_system()

    # Parameters
    e = [0, 1]        # Jones polarisation vector
    θ, ϕ = 0.0, 0.0   # rotation angles [rad]

    k = 2π / lambda
    ε = model(lambda)

    kvec = [0, 0, k]
    E₀ = [e..., 0]
    alpha = LDR(ε, grid.spacing, kvec, E₀)
    alphas = fill(alpha, size(inds))

    # 3. calculated the incident field Einc, at each dipole,
    pw = PlaneWave(k, e, θ, ϕ)
    Eincs = [field(pw, coord) for coord in coords]

    # 5. Define the DDA problem
    prob = GridProblem(k, norm(e), grid, inds, alphas, Eincs)

    # 6. Solve the DDA problem
    sol = solve(prob, BiCGStabl(), reltol = 1e-4, verbose = verbose)
    return C_abs(sol), C_sca(sol)

end
#
#
# to = TimerOutput()
# @timeit to "main" begin
#     Threads.@threads :static for i = 1:length(lambda_range)
#         @show i
#
#         Q_abs[i], Q_sca[i] = run_calc(lambda_range[i], model)
#
#     end
# end
# print_timer(to)

to = TimerOutput()
for i = eachindex(lambda_range)
    @timeit to "main" begin
        @show i

        Q_abs[i], Q_sca[i] = run_calc(lambda_range[i], model)
    end
    print_timer(to)
end



layout = Layout()
trace = [
    scatter(x=lambda_range, y=Q_abs, name = "abs"),
    scatter(x=lambda_range, y=Q_sca, name = "sca"),

]
plot(trace, layout)



plot(lambda_range, Q_abs, label = "abs");
plot!(lambda_range, Q_sca, label = "sca");
plot!(lambda_range, Q_abs + Q_sca, label = "ext")




lambda_range = 620:1.0:740


lamda_range









# target: Au disk
# d = 10 nm
# # of dipoles along the diameter: 50
# heightt = 1 nm (5 dipoles)
# grid spacing: d = 0.2 nm
# # of dipoles: ??
# gap: d = (1) 0.2 nm, (2) 0.4 nm, (5) 1 nm,  (100) 20 nm
# >> peak resonance at ??? nm
# # of dipoles: 9880 ?
# wavelength range: 620 740



# 1. Define a grid
# disk_r = 5. # nm


diameter = 10.0 # nm
height = 1.0 # nm

gap = 1.0 # nm
spacing = 0.2
# sum(occ) = 19760

gap = 4 / 3 # nm
spacing = 1 / 3
# sum(occ) = 4296

Nx = round(Int, (2 * diameter + gap) / spacing)
Ny = round(Int, diameter / spacing)
Nz = round(Int, height / spacing)

offset = [spacing / 2, spacing / 2, spacing / 2]
origin = [-(2 * diameter + gap) / 2, -diameter / 2, -height / 2] .+ offset

grid = CartesianGrid(origin, [spacing, spacing, spacing], (Nx, Ny, Nz))
grid

# 2. Define the target(s)
origin = [(diameter + gap) / 2, 0, 0]

d1 = DDA.Disk(-origin, diameter / 2, height)
d2 = DDA.Disk(origin, diameter / 2, height)
s = DDA.Composite([d1, d2])

occ = DDA.discretize(grid, s)

layout = Layout(yaxis = attr(scaleanchor = "x", scaleratio = 1))
plot(heatmap(z = occ[:, :, 1]' .* 1, colorscale = "Viridis", xgap = 3, ygap = 3), layout)

coords = grid[occ]

trace = scatter3d(x = [p[1] for p in coords], y = [p[2] for p in coords], z = [p[3] for p in coords], marker = attr(color = 1:length(coords)), mode = "markers")
layout = Layout(scene = attr(aspectmode = :data, camera_projection_type = :orthographic))
plot(trace, layout)


# 3. Define material

# All in one together


using Unitful
using PhysicalConstants.CODATA2018: ε_0, c_0, h

using RefractiveIndexDatabase
using Permittivity

λ_range = (632.0:1.0:732.0)u"nm"
Au = get_material("main", "Au", "Babar")

permittivity(material) = @. (material.n + material.k * im)^2
# permittivity(material) = @. (material.n^2 - material.k^2) + (2 * material.n * material.k) * im
wavelength(material) = material.λ * u"μm"

ε = permittivity(Au)
λ = wavelength(Au) .|> u"nm"

int = QuadraticInterpolation(ε, λ)
ε_int = int.(λ_range)

#trace1 = scatter(x=ω_range, y=real(ε_int), name="Re(ε(ω))")
trace2 = scatter(x = ω_range, y = imag(ε_int), name = "Im(ε(ω))")
trace3 = scatter(x = ω, y = real(ε), name = "Re(ε)", mode = "markers")
trace4 = scatter(x = ω, y = imag(ε), name = "Im(ε)", mode = "markers")

layout1 = Layout(xaxis = attr(title = "[eV]", range = extrema(ω_range)), yaxis = attr(title = "ε", range = [-20, 0]))
layout2 = Layout(xaxis = attr(title = "[eV]", range = extrema(ω_range)), yaxis = attr(title = "ε", range = [0, 1]))

[plot([trace1, trace3], layout1) plot([trace2, trace4], layout2)]





# 4. Define incindent field
k = 2π     # wavenumber
e = [1, 0]  # Jones polarisation vector
θ, ϕ = 0.0, 0.0 # rotation angles [rad]


E_inc = DDA.PlaneWave(k, e, θ, ϕ)
E = DDA.field(E_inc, grid[occ])



