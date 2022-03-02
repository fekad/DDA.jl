using DDA

# using Unitful
# using PhysicalConstants.CODATA2018: c_0
# using RefractiveIndexDatabase
# using Permittivity

using StaticArrays
using IterativeSolvers
using LinearAlgebra

using BenchmarkTools

using Plots
plotlyjs()


############################################################################
# High-level interface for DDA calculations
############################################################################
# grid size => number of dipoles
# 8 => 304
# 12 => 1064
# 16 => 2320
# 24 => 7664
# 32 => 17904
# 48 => 59728

# k = 2π/λ
# |m|kd = 1
# x = ka =2πa/λ

k = 2π
kvec = [0, 0, k] # wavevector
E₀  = ComplexF64[1, 0, 0]   # polarisation vector

ε = 1.33 + 0.1im

# 1. create the coordinates of the dipoles

origin = [0.0, 0.0, 0.0]
spacing = 1.0
dims = (8, 8, 8)

g = CartesianGrid(origin, spacing, dims)

origin = DDA.center(g)
radius = (8 + 0.49) / 2

sphere = Sphere(origin, radius)
coords = coordinates(g, sphere)

# 2. assign the polarizability αj to each dipole,

alphas = fill(LDR(ε, spacing, kvec, E₀), size(coords))

# 3. calculated the incident field Einc,j at each dipole,

pw = PlaneWave(SVector(kvec...), SVector(E₀...))
Eincs = [field(pw, coord) for coord in coords ]

# 4. assemble the interaction matrix A and
prob = DipoleProblem(k, norm(E₀), coords, alphas, Eincs)

A = DDA.interactions(k, coords, alphas)

# 5. solve for P in the system of linear equations
# P = minres(A, reinterpret(ComplexF64,Eincs), verbose=:true)
# P = bicgstabl(A, reinterpret(ComplexF64,Eincs), verbose=:true)
sol = solve(prob, BiCGStabl(), reltol = 1e-8)


@show C_abs(sol)
@show C_sca(sol)
@show C_ext(sol)
# C_abs(sol) = 1.2892547488225787

############################################################################
# High-level interface for DDA calculations
############################################################################


function sphere_system(Nd, k, m)

    # Parameters
    d = 1.0           # spacing
    e = [1, 0]        # Jones polarisation vector
    θ, ϕ = 0.0, 0.0   # rotation angles [rad]

    grid = CartesianGrid([0.0, 0.0, 0.0], d, (Nd, Nd, Nd))
    sphere = Sphere(DDA.center(grid), (Nd + 0.49) / 2)

    # 1. create the coordinates of the dipoles,
    coords = coordinates(g, sphere)

    # 2. assign the polarizability αj to each dipole
    # scatterer = Scatterer(sphere, LDRModel(m^2))
    alphas = fill(LDR(m^2, d, kvec, E₀), size(coords))

    # 3. calculated the incident field Einc, at each dipole,
    pw = PlaneWave(k, e, θ, ϕ)
    Eincs = [field(pw, coord) for coord in coords ]

    # 5. Define the DDA problem
    prob = DipoleProblem(k, norm(E₀), coords, alphas, Eincs)

    # 6. Solve the DDA problem
    return solve(prob, BiCGStabl(), reltol = 1e-8)

end

# Parameters:
Nd = 8
a = Nd / 2

m = 1.33 + 0.01im

k = range(0., 12.5 / a, length = 101)[2:end]
Q_abs = zeros(length(k), 1)
Q_sca = zeros(length(k), 1)

for (i, k) = enumerate(k)
    local sol
    sol = sphere_system(Nd, k, m)
    Q_abs[i] = C_abs(sol) / (π * a^2)
    Q_sca[i] = C_sca(sol) / (π * a^2)
end

plot!(k * a, Q_abs, label = "abs", yscale = :log10, ylim = [0.005, 5])
plot!(k * a, Q_sca, label = "sca")


