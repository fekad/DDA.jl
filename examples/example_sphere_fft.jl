using DDA


# using Unitful
# using PhysicalConstants.CODATA2018: c_0
# using RefractiveIndexDatabase
# using Permittivity

using StaticArrays
using IterativeSolvers
using LinearAlgebra: norm

using BenchmarkTools

using Plots
plotlyjs()


############################################################################
# High-level interface for DDA calculations
############################################################################


k = 2π      # wavenumber
e = [1, 0]  # Jones polarisation vector
θ, ϕ = 0.0, 0.0 # rotation angles [rad]

ε = 1.33 + 0.1im

# 1. Define a grid
origin = [0.0, 0.0, 0.0]
spacing = 1.0
dims = (32, 32, 32)

grid = CartesianGrid(origin, spacing, dims)

# 2. Define the target(s)
origin = DDA.center(grid)
radius = 5.0

sphere = Sphere(origin, radius)

# coords, occ, alphas = discretize(p)
occ = discretize(grid, sphere)
inds = indices(grid, sphere)
coords = grid[inds]

# 3. Define the material properties
# model = LDRModel(ε)
# scatterer = Scatterer(sphere, model)

pw = PlaneWave(k, e, θ, ϕ)
alphas = fill(LDR(ε, spacing, pw.kvec, pw.E₀), size(coords))


# 4. Define incindent field
# pw = PlaneWave(SVector(kvec...), SVector(E₀...))
# Eincs = similar(coords, SVector{3,Complex{Float64}})
# for i = eachindex(coords)
#     Eincs[i] = field(p.Einc, coords[i])
# end
Eincs = [field(pw, coord) for coord in coords]


# 5. Define the DDA problem
# prob = GridProblem(grid, scatterer, Einc)
prob = GridProblem(k, norm(e), grid, inds, alphas, Eincs)


# 6. Solve the DDA problem
sol = solve(prob, BiCGStabl(), reltol = 1e-5)


@show C_abs(sol)
@show C_sca(sol)
@show C_ext(sol)


############################################################################
# High-level interface for DDA calculations
############################################################################

function sphere_system(Nd, k, m)

    # Parameters
    d = 1.0           # spacing
    e = [1, 0]        # Jones polarisation vector
    θ, ϕ = 0.0, 0.0   # rotation angles [rad]

    # 1. create the coordinates of the dipoles,
    grid = CartesianGrid([0.0, 0.0, 0.0], d, (Nd, Nd, Nd))
    sphere = Sphere(DDA.center(grid), (Nd + 0.49) / 2)

    # 2. assign the polarizability αj to each dipole
    # scatterer = Scatterer(sphere, LDRModel(m^2))

    # occ = discretize(grid, sphere)
    # coords = grid[occ]
    inds = indices(grid, sphere)
    coords = grid[inds]

    ε = m^2
    kvec = [0, 0, k]
    E₀ = [e..., 0]
    alphas = fill(LDR(ε, spacing, kvec, E₀), size(coords))

    # 3. calculated the incident field Einc, at each dipole,
    # Einc = PlaneWave(k, e, θ, ϕ)
    pw = PlaneWave(k, e, θ, ϕ)
    Eincs = [field(pw, coord) for coord in coords]

    # 5. Define the DDA problem
    prob = GridProblem(k, norm(e), grid, inds, alphas, Eincs)

    # 6. Solve the DDA problem
    return solve(prob, BiCGStabl(), reltol = 1e-8)

end

# Parameters:
Nd = 24
Nd = 8
a = Nd / 2

m = 1.33 + 0.01im

k = range(0.0, 12.5 / a, length = 101)[2:end]
Q_abs = zeros(length(k), 3)
Q_sca = zeros(length(k), 3)

for (i, k) = enumerate(k)
    local sol
    sol = sphere_system(Nd, k, m)
    Q_abs[i] = C_abs(sol) / (π * a^2)
    Q_sca[i] = C_sca(sol) / (π * a^2)
end

plot(k * a, Q_abs, label = "abs", yscale = :log10, ylim = [0.005, 5]);
plot!(k * a, Q_sca, label = "sca")



############################################################################
# "Mid-level" interface to DDA calulations
############################################################################

function sphere_system_steps(Nd, k, m)

    # Parameters
    d = 1.0          # spacing
    e = [1, 0]      # Jones polarisation vector
    θ, ϕ = 0.0, 0.0   # rotation angles [rad]

    # 1. create the coordinates of the dipoles,
    grid = CartesianGrid([0.0, 0.0, 0.0], d, (Nd, Nd, Nd))

    # 2. Define and discretise the target(s)
    origin = DDA.center(grid)
    radius = (Nd + 0.49) / 2

    t = Sphere(origin, radius)

    inds = indices(grid, t)
    occ = discretize(grid, t)
    coords = grid[occ]

    N = length(coords)
    a_eff = (N * d^3 * 3 / 4π)^(1 / 3)
    @show length(coords)
    @show radius / a_eff * 100

    # 3. Assign the polarizability to and calculate the incident field for each dipole
    # |m|kd
    @show abs(m) * k * d

    ε = m^2
    kvec = [0, 0, k]
    E₀ = [e..., 0]

    α = LDR(ε, d, kvec, E₀)
    alphas = fill(α, length(coords))

    pw = PlaneWave(k, e, θ, ϕ)

    Einc = similar(coords, SVector{3,Complex{Float64}})
    for i = eachindex(coords)
        Einc[i] = field(pw, coords[i])
    end

    # 4. Assemble the interaction matrix A
    # The solution using FFT
    # A = InteractionTensor(g, k, α, occ)
    # P = bicgstabl(factorise(A), E; abstol=1e-8)
    A = DDA.TensorConvolution(grid, inds, k, alphas)

    # memory requirements:
    # (3*N)^2 * 2 * 8 / 1024^3

    # 5. solve the system of linear equations to get the polarisabilities
    P = similar(coords, SVector{3,Complex{Float64}})
    fill!(P, zero(SVector{3,Complex{Float64}}))

    bicgstabl!(reinterpret(ComplexF64, P), A, reinterpret(ComplexF64, Einc); reltol = 1e-4, verbose = true)

    # 6. Analyse the results

    return C_abs(k, E₀, P, alphas), C_ext(k, E₀, Einc, P), C_sca(k, E₀, Einc, P, alphas)

end


# Pararmeters:
Nd = 24
a = Nd / 2

# m = 0.96 + 1.01im
m = 1.33 + 0.01im
# m = 2 + im

k = range(0.0, 12.5 / a, length = 101)[2:end]

out = zeros(length(k), 3)
for i = eachindex(k)
    out[i, :] .= sphere_system_steps(Nd, k[i], m)
    @show out[i, :]
end

plot(k * a, out[:, 1] / (π * a^2), label = "abs", yscale = :log10, ylim = [0.005, 5]);
plot!(k * a, out[:, 3] / (π * a^2), label = "sca")


