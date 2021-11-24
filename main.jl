using Revise
using DDA



############################################################################
# High-level interface for DDA calcularions
############################################################################

# 1. Define a grid

Nx, Ny, Nz = 10, 10, 10
origin = [0, 0, 0]
spacing = [1, 1, 1]

g = CartesianGrid(origin, spacing, (Nx, Ny, Nz))
# g = DDA.CartesianGrid(origin, spacing, [Nx, Ny, Nz]) fails!


# 2. Define the target(s)

center = [0,0,0]
radius = 1

t = Sphere(center, radius)
# t = Sphere(grid, center, radius) ???

# 3. Define the material properties
DDA.MaterialModels.test
# ε = Material.gold("") ???
ε = Material.load("Ag.yaml")
α = LDR(ε)

s = Scatterer(t, α)


# 4. Define the DDA problem

k = 2π      # wavenumber
e = [1, 0]  # Jones polarisation vector
θ, ϕ = 0, 0 # rotation angles [rad]

E_inc = PlaneWave(k, e, θ, ϕ)

a = 1 # effective radii ??

prob = DDAGridProblem(grid, a, s, E_inc)


# 5. Solve the DDa problem

sol = solve(prob, BiCGStablFFT(), tol=1e-12)


############################################################################
# Notes:
############################################################################

# Analyse the results

@show C_abs(sol)
@show C_ext(sol)
@show C_sca(sol)

# calculate the field at a point
r = [1, 1, 1]
E = E_sca_FF(sol, r)

# show the discretised target
plot(grid, t)

# real and imaginary part of permittivity
plot(ε, LinRange(10, 100, 1001))

# E field in 3D
plot(sol)

# get the taylor series at a point

r = [1, 1, 1]
order = 3
s = taylor(sol, r, order)

# use case: calulation at different wavelengths

krange = LinRange(100, 500, 1001)
Q_abs = zeros(length(krange))

for (i, k) in enumerate(krange)

    E_inc = PlaneWave(k, e, 0, 0)
    prob = DDAGridProblem(grid, (t, α), E_inc)

    sol = solve(prob, BiCGStablFFT(), tol=1e-12)

    filename = "dda_calc_$i.sav"
    save(sol, filename)

    Q_abs[i] = C_abs(sol)

end

plot(krange, Q_abs)


# dipoles without grid

# dipoles [Nx3]
# polarisabilities [N] or [Nx3x3]
# E_inc [Nx3]
prob = DDAProblem(dipoles, polarisabilities, E_inc)


# polarizability_LDR(eps, d, k, E₀) ???
# grid: d, E_inc: k, E₀


############################################################################
# "Medium"" level interface to DDA calulations
############################################################################

# everything stored as a list and there is a mapping between idicies
# r (Nd, 3)
# E (Nd, 3)
# α (Nd, 3)
# inds (Nd, 3)
# A (Nx, Ny, Nz, 6) ??? lazy?


# parameters
k = 2π      # wavenumber
a = 1       # effective radii
e = [1, 0]  # Jones polarisation vector
θ, ϕ = 0, 0 # rotation angles [rad]


# steps:
# 1. create the coordinates of the dipoles,

Nx, Ny, Nz = 10, 10, 10
origin = [0, 0, 0]
spacing =  [1, 1, 1]

g = Grid(origin, spacing, [Nx, Ny, Nz])

center = [0,0,0]
radius = 1

t = Sphere(center, radius)
r, inds = dipoles(g, s)

# a = effective radii ???

# 2. assign the polarizability αj to each dipole

ε = permitivity(Material.load("Ag.yaml"), k)
d = a * spacing # ???
E₀ = norm(e)
α = LDR(ε, d, k, E₀)


# 3. calculated the incident field E_inc, at each dipole,

E_inc = field(PlaneWave(k, e, θ, ϕ), r)


# 4. assemble the interaction matrix A and

A = interactions(k, r, α)

# 5. solve for P in the system of linear equations

P = A \ E
P = bicgstabl(A, E; abstol=1e-8)
P = cg(A, E; abstol=1e-8, maxiter=1000, verbose=true)
P = minres(A, E; abstol=1e-8, maxiter=1000)


# The solution using FFT
A = InteractionTensor(g, k, α, inds)
P = bicgstabl(factorise(A), E; abstol=1e-8)


# Analyse the results

@show C_abs(k, P, α)
@show C_ext(k, E0, Ei, P)
@show C_sca(k, E0, Ei, P, α)

# calculate the field at a point
r = [1, 1, 1]
E = E_sca_FF(k, r, P, r_E)


end

# using Plots
#
# plotlyjs()
#
# # k = [0, 0, 1.]
# #
# #
# # E_0 = 1.
# #
# # E_inc = DDA.PlaneWave(E_0, k)
# #
# #
# # r_n = [1, 1, 1]
# # r = r_n * (0:0.1:10)'
# #
# # E = DDA.field(r[:,1], E_inc)
# # E = DDA.field(r, E_inc)
# #
# #
# # plot(real(E))
#
#
# # target: Au disk
# # d = 10 nm
# # # of dipoles along the diameter: 50
# # grid spacing: d = 0.2 nm
# # # of dipoles: ??
# # spaceing: d = (1) 0.2 nm, (2) 0.4 nm, (5) 1 nm,  (100) 20 nm
# # >> peak resonance at ??? nm
#
# # steps:
# # 1. create the coordinates of the dipoles,
# # 2. assign the polarizability αj to each dipole,
# # 3. calculated the incident field E_inc, at each dipole,
# # 4. assemble the interaction matrix A and
# # 5. solve for P in the system of linear equations
#
# # NOTE: the dimension of the grid shouldn't matter.
# #  - but it hard to wald throught the whole space
# #  - based on the scatterer create the grid using extremas (like bounding box)
# g = CartesianGrid([-20.5,-20.5,-20.5], ones(3), (41, 41, 41))
# for i in 1:16
#     s = Sphere(i, [0, 0, 0])
#     r = dipoles(g, s)
#     @show i, length(r)
# end
#
#
#
# # 1. create the coordinates of the dipoles,
#
# g = CartesianGrid([1/6,1/6,1/6], 1/3*ones(3), (30, 30, 30))
# s = Sphere(5., [5, 5, 5])
# r = dipoles(g, s)
#
# rr = reduce(hcat, r)
# scatter3d(rr[1,:],rr[2,:],rr[3,:])
#
#
#
# g = CartesianGrid([1.5,1.5,1.5], 3*ones(3), (30, 30, 30))
# s = Sphere(15., [15, 15, 15])
# r = dipoles(g, s)
#
# rr = reduce(hcat, r)
# scatter3d(rr[1,:],rr[2,:],rr[3,:])
#
# # g = CartesianGrid([-14.5, -14.5, -14.5], ones(3), (30,30,30))
# # s = Sphere(15., [0, 0, 0])
# # d = dipoles(g,s)
#
#
#
# g = CartesianGrid([-100, -100, -100], 200*ones(3), (11, 11, 11))
# s = Disk(15, 5, [15, 15, 0 ])
# r = dipoles(g, s)
# rr = reduce(hcat, r)
#
#
#
# g = CartesianGrid([-200, -200, -200], 40*ones(3), (11, 11, 11))
# s = Sphere(200., [0, 0, 0])
# r = dipoles(g, s)
# rr = reduce(hcat, r)
# # gr()
# plotlyjs()
# scatter(rr[1,:],rr[2,:],rr[3,:], aspect_ratio=:auto)
# xlims!(-200, 200); ylims!(-200, 200); zlims!(-200, 200)
#
#
#
#
#
# # target: Au disk
# # d = 10 nm
# # height = 1 nm (5 dipoles)
# # # of dipoles along the diameter: 50
# # grid spacing: d = 0.2 nm
# # spaceing: d = (1) 0.2 nm, (2) 0.4 nm, (5) 1 nm,  (100) 20 nm
# # >> peak resonance at ??? nm
# # # of dipoles: 9880?
#
#
# d = 1/3
# g = CartesianGrid([d/2, d/2, 0], d*ones(3), (30, 30, 3))
# s = Disk(5., 1, [5, 5, 0])
# r = dipoles(g, s)
# rr = reduce(hcat, r)
#
# plotlyjs()
# scatter(rr[1,:],rr[2,:],rr[3,:], aspect_ratio=:auto)
# xlims!(0, 10); ylims!(0, 10); zlims!(0, 10)
#
#
# g = CartesianGrid([0.1, 0.1, 0], [.2, .2, .2], (50, 50, 5))
# s = Disk(5, 1, [5, 5, 0 ])
# r = dipoles(g, s)
# rr = reduce(hcat, r)
#
# # gr()
# plotlyjs()
# scatter(rr[1,:],rr[2,:],rr[3,:], aspect_ratio=:auto)
# xlims!(0, 10); ylims!(0, 10); zlims!(0, 10)
#
#
#
#
#
#
# # 2. assign the polarizability αj to each dipole,
#
# E₀ = [1, 0, 0]
# k = 2π * [0, 0, 1]
#
# n = 4.686878048780487+0.08797560975609756im
# eps = n^2
# d = g.spacing[1]
#
#
# alpha = polarizability_LDR(eps, d, k, E₀)
#
#
# # 3. calculated the incident field E_inc, at each dipole,
#
# p = PlaneWave(E₀, k)
# E_inc = field(r, p)
#
#
# # 4. assemble the interaction matrix A and
# using LinearAlgebra
#
# @time A = DDA.interaction(norm(k), r, alpha * ones(length(r)))
#
#
#
# # 5. solve for P in the system of linear equations
#
#
# using IterativeSolvers
#
# @time s1=cg(A, E_inc[:]; abstol=1e-8, maxiter=1000, verbose=true)
# @time s2=minres(A, E_inc[:]; abstol=1e-8, maxiter=10000)
# @time s3=bicgstabl(A, E_inc[:]; abstol=1e-8)
#
# @time s4=A \ E_inc[:]