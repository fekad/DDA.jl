using DDA

# using Unitful
# using PhysicalConstants.CODATA2018: c_0
# using RefractiveIndexDatabase
# using Permittivity

using StaticArrays
using IterativeSolvers

using BenchmarkTools

using Plots
plotlyjs()


############################################################################
# High-level interface for DDA calculations
############################################################################

# 1. Define a grid
origin = [0.0, 0.0, 0.0]
spacing = [1.0, 1.0, 1.0]
dims = (32, 32, 32)

grid = CartesianGrid(origin, spacing, dims)

# 2. Define the target(s)
origin = DDA.center(grid)
radius = 5.0

sphere = DDA.Sphere(origin, radius)

# 3. Define the material properties
ε = 1.33 + 0.1im
model = DDA.LDRModel(ε)

scatterer = DDA.Scatterer(sphere, model)

# 4. Define incindent field
k = 2π      # wavenumber
e = [1, 0]  # Jones polarisation vector
θ, ϕ = 0.0, 0.0 # rotation angles [rad]

Einc = PlaneWave(k, e, θ, ϕ)

# 5. Define the DDA problem
prob = DDA.GridProblem(grid, scatterer, Einc)

# 6. Solve the DDA problem
sol = DDA.solve(prob, DDA.BiCGStablFFT(), tol = 1e-5)


@show DDA.C_abs(sol)


############################################################################
# High-level interface for DDA calculations
############################################################################

function sphere_system(Nd, k, m)

    # Parameters
    d = 1.0           # spacing
    e = [1, 0]        # Jones polarisation vector
    θ, ϕ = 0.0, 0.0   # rotation angles [rad]

    # 1. create the coordinates of the dipoles,
    grid = CartesianGrid([0.0, 0.0, 0.0], [d, d, d], (Nd, Nd, Nd))
    sphere = DDA.Sphere(DDA.center(grid), (Nd + 0.49) / 2)

    # 2. assign the polarizability αj to each dipole
    scatterer = DDA.Scatterer(sphere, DDA.LDRModel(m^2))

    # 3. calculated the incident field Einc, at each dipole,
    Einc = DDA.PlaneWave(k, e, θ, ϕ)

    # 5. Define the DDA problem
    prob = DDA.GridProblem(grid, scatterer, Einc)

    # 6. Solve the DDA problem
    return DDA.solve(prob, DDA.BiCGStablFFT(), tol = 1e-12)

end

# Pararmeters:
Nd = 24
a = Nd / 2

m = 1.33 + 0.01im

k = range(0., 12.5 / a, length = 101)[2:end]
Q_abs = zeros(length(k), 3)
# Q_sca = zeros(length(k), 3)

for (i, k) = enumerate(k)
    local sol
    sol = sphere_system(Nd, k, m)
    Q_abs[i] = DDA.C_abs(sol) / (π * a^2)
    # Q_sca[i] = DDA.C_sca(sol) / (π * a^2)
end

plot(k * a, Q_abs, label = "abs", yscale = :log10, ylim = [0.005, 5]);
# plot!(k * a, Q_sca, label = "sca")



############################################################################
# "Mid-level" interface to DDA calulations
############################################################################

function sphere_system_steps(Nd, k, m)

    # Parameters
    d = 1.0          # spacing
    e = [1, 0]      # Jones polarisation vector
    θ, ϕ = 0.0, 0.0   # rotation angles [rad]

    # 1. create the coordinates of the dipoles,
    g = CartesianGrid([0.0, 0.0, 0.0], [d, d, d], (Nd, Nd, Nd))

    # 2. Define and discretise the target(s)
    origin = DDA.center(g)
    radius = (Nd + 0.49) / 2

    t = DDA.Sphere(origin, radius)

    occ = DDA.discretize(g, t)
    coords = g[occ]

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

    pw = DDA.PlaneWave(k, e, θ, ϕ)

    Einc = similar(coords, SVector{3,Complex{Float64}})
    for i = eachindex(coords)
        Einc[i] = field(pw, coords[i])
    end

    # 4. Assemble the interaction matrix A
    # The solution using FFT
    # A = InteractionTensor(g, k, α, occ)
    # P = bicgstabl(factorise(A), E; abstol=1e-8)
    A_conv = TensorConvolution(g, occ, k, α)

    # memory requirements:
    # (3*N)^2 * 2 * 8 / 1024^3

    # 5. solve the system of linear equations to get the polarisabilities
    P = similar(coords, SVector{3,Complex{Float64}})
    fill!(P, zero(SVector{3,Complex{Float64}}))

    bicgstabl!(reinterpret(ComplexF64, P), A_conv, reinterpret(ComplexF64, Einc); reltol = 1e-4, verbose = true)

    # 6. Analyse the results

    return DDA.C_abs(k, E₀, P, alphas[1]), DDA.C_ext(k, E₀, Einc, P), DDA.C_sca(k, E₀, Einc, P, alphas[1])

end


# Pararmeters:
Nd = 24
a = Nd / 2

# m = 0.96 + 1.01im
m = 1.33 + 0.01im
# m = 2 + im

k = range(0., 12.5 / a, length = 101)[2:end]

out = zeros(length(k), 3)
for i = eachindex(k)
    out[i, :] .= sphere_system_steps(Nd, k[i], m)
    @show out[i, :]
end

plot(k * a, out[:, 1] / (π * a^2), label = "abs", yscale = :log10, ylim = [0.005, 5]);
plot!(k * a, out[:, 3] / (π * a^2), label = "sca")



# TODO:
# - Multiple scatterers:
#
#     # 3. Define the material properties
#     sphere = DDA.Sphere(DDA.center(grid), 5.)
#     single_dipole = DDA.Dipole([1., 1., 1.])
#
#     ε1 = 1.33 + 0.1im
#     model1 = DDA.LDRModel(ε1)
#     scatterer1 = DDA.Scatterer(sphere, model1)
#
#     ε2 = Lorentz(params...)(k)
#     model2 = Direct(ε2) # DDA.Atomic(ε2)
#     scatterer2 = DDA.Scatterer(single_dipole, model2)
#
#     scatterers = [scatterer1, scatterer2]
#
# - Define the material properties as function/model
#     # TODO: wavelelngth or frequency???
#     # TODO: refractive index or permittivity???
#
#     Ag = get_material("main","Ag","Johnson")
#
#     λ_Ag = Ag.λ
#     ñ_Ag = Ag.n + Ag.k * im
#     ε_Ag = ñ_Ag.^2;
#     ω_Ag = ustrip.(2π * c_0 ./ (λ_Ag * 1u"μm") .|> u"THz")
#
#     ε = PermittivityTable(ω_Ag, ε_Ag)
#
#     # ε = Material.gold("") ???
#     # ε = Material.load("Ag.yaml")
#     #
#     α = LDR(ε)
#     s = Scatterer(t, α)
#
# - Random notes:
#     # grid size => number of dipoles
#     # 8 => 304
#     # 12 => 1064
#     # 16 => 2320
#     # 24 => 7664
#     # 32 => 17904
#     # 48 => 59728
#
#     # k = 2π/λ
#     # |m|kd = 1
#     # x = ka =2πa/λ
#
#
#     # 5. solve for P in the system of linear equations
#
#     Ei = reinterpret(ComplexF64, Einc)
#     Pi = A \ Ei
#     Pi = bicgstabl(A, Ei; abstol = 1e-3, verbose = true)
#     # Pi = cg(A, Ei; abstol=1e-3, maxiter=1000, verbose=true)
#     # Pi = minres(A, Ei; abstol=1e-3, maxiter=1000, verbose=true)
#
#     P = reinterpret(SVector{3,ComplexF64}, Pi)
#
#
#     # Analyse the results
#
#     @show C_abs(sol)
#     @show C_ext(sol)
#     @show C_sca(sol)
#
#     # calculate the field at a point
#     r = [1, 1, 1]
#     E = E_sca_FF(sol, r)
#
#     # show the discretised target
#     plot(grid, t)
#
#     # real and imaginary part of permittivity
#     plot(ε, LinRange(10, 100, 1001))
#
#     # E field in 3D
#     plot(sol)
#
#     # get the taylor series at a point
#
#     r = [1, 1, 1]
#     order = 3
#     s = taylor(sol, r, order)
#
#
#     # dipoles without grid
#
#     # dipoles [Nx3]
#     # polarisabilities [N] or [Nx3x3]
#     # Einc [Nx3]
#     prob = DDAProblem(dipoles, polarisabilities, Einc)
#
