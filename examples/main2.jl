# using Revise
using DDA
using RefractiveIndexDatabase
using Permittivity
using Unitful
using BenchmarkTools
using PhysicalConstants.CODATA2018: ε_0, c_0

using StaticArrays
using IterativeSolvers

using Plots
plotlyjs()

# DDSCAT
# The units must be the same as the wavelength units used in the file specifying the refractive index.
# (x, y, z)_n = [(i, j, k)_n     + (x0, y0, z0)] × d, where d is the lattice constant (in physical units)
# The dipole spacing d in physical units is determined from the specified value of aeff and the number N of dipoles in the target: d = (4π/3N)1/3aeff.






# Notes for ELLIPSOID example.
#
# ddscat.par is set up to calculate absorption and scattering by
# a ellipsoid with diameters/d = 48.49 48.49 48.49  (d=interdipole separation).
# this of course means a sphere of diameter = 48.49*d
#
# m = 0.96+1.01i (refractive index of Au at 500 nm)
# effective radius = 5 units
# lambda = 6.283185 units  (wavelength in vacuum)
# this gives
# size parameter x = 2*pi*a/lambda = 2*pi*5/2*pi = 5
#
# TARGET --- Sphe       e,  59728 dipoles, 1.000 1.000 1.000=x,y,z lattice spacing
# 0.04123856 = d/aeff for this target [d=dipole spacing]
#   0.016408 = d (physical units)
# ----- physical extent of target volume in Target Frame ------
#  -0.393802      0.393802 = xmin,xmax (physical units)
#  -0.393802      0.393802 = ymin,ymax (physical units)
#  -0.393802      0.393802 = zmin,zmax (physical units)
# AEFF=      0.397890 = effective radius (physical units)
# WAVE=      0.500000 = wavelength (in vacuo, physical units)
# K*AEFF=      5.000033 = 2*pi*aeff/lambda


# k = 2π / λ /10 # wavenumber
k=1.
# function runtest(k)

    # parameters
    λ = .5      # wavelength (microns)
    k = 2π / λ  # wavenumber
    a = 0.39789 # effective radii (microns)
    d = 1.      # spacing
    Nd = 32
    Nd = 24
    e = [1, 0]  # Jones polarisation vector
    θ, ϕ = 0, 0 # rotation angles [rad]

    # steps:
    # 1. create the coordinates of the dipoles,
    # CartesianGrid{Float64,3}(SVector(0.,0.,0.),SVector(1.,1.,1.),Dims((11,11,11)))
    # CartesianGrid{Float64,3}((0.,0.,0.),(1.,1.,1.),(11,11,11))

    Nx = Ny = Nz = Nd
    origin = Point3(0., 0., 0.)

    g = CartesianGrid(origin, [d, d, d], (Nx, Ny, Nz))

    center = DDA.center(g)
    radius = (Nd + .49)/2

    t = DDA.Sphere(center, radius)

    r2, occ = DDA.dipoles(g, t)
    occ = DDA.discretize(g, t)
    # r2 == g[occ]

    N = length(r2)

    # a = effective radii ???
    a_eff = (3*N*(a*d)^3/(4*pi))^(1/3)



    a * d / a_eff       # d/aeff for this target [d=dipole spacing]
    a * d / a_eff * a   # d (physical units) ???

    # ----- physical extent of target volume in Target Frame ------
    #      -0.393802      0.393802 = xmin,xmax (physical units)
    #      -0.393802      0.393802 = ymin,ymax (physical units)
    #      -0.393802      0.393802 = zmin,zmax (physical units)
    Nd/2 * (a * d / a_eff * a)

    a               # effective radius (physical units)
    λ               # wavelength (in vacuo, physical units)
    k*a             # 2*pi*aeff/lambda


    # ( 0.20619  0.00000  0.00000 ) = k vector (latt. units) in TF
    k* (a * d / a_eff * a)





    # 2. assign the polarizability αj to each dipole

    # m = 0.96 + 1.01im
    m = 1.33 + 0.01im

    ε = m^2

    # |m|kd
    abs(m) * k * (a * d / a_eff * a)


    # ε = permitivity(Material.load("Ag.yaml"), k)

    # d = a * spacing # ???
    kvec = [0, 0, k]
    E₀ = [e..., 0]
    α = LDR(ε, d, kvec, E₀)

    alphas = fill(α, length(r2))


    # 3. calculated the incident field E_inc, at each dipole,
    pw = DDA.PlaneWave2(k, e, θ, ϕ)
    E_inc = field(pw, r2)



    # 4. assemble the interaction matrix A and

    # memory requirements:
    (3*N)^2 * 2 * 8 / 1024^3


    # 5. solve for P in the system of linear equations


    A_conv = TensorConvolution(g, occ, k, α)

    Ei = reinterpret(ComplexF64, E_inc)
    P3 = bicgstabl(A_conv, Ei; abstol=1e-4, verbose=true)

    P = reinterpret(SVector{3,ComplexF64}, P3)





    # # The solution using FFT
    # A = InteractionTensor(g, k, α, occ)
    # P = bicgstabl(factorise(A), E; abstol=1e-8)


    # Analyse the results

    # @show DDA.C_abs(k, E₀, P, alphas)
    # @show DDA.C_ext(k, E₀, E_inc, P)
    # @show DDA.C_sca(k, E₀, E_inc, P, alphas)
    @show k
    return DDA.C_abs(k, E₀, P, alphas), DDA.C_ext(k, E₀, E_inc, P), DDA.C_sca(k, E₀, E_inc, P, alphas)
end
# return DDA.C_abs(k, E₀, P, alphas), DDA.C_ext(k, E₀, E_inc, P), DDA.C_sca(k, E₀, E_inc, P, alphas)



k = range(0,1,length=101)[2:end]
out = zeros(length(k), 3)
for i = 1:length(k)
    out[i,:] .= runtest(k[i])
end



plot(k*13, out[:,1], label="abs", yscale=:log10, ylim=[1e0, 1.7e3])
plot!(k*13, out[:,3], label="sca")

nothing































# grid size => number of dipoles
# 8 => 304
# 12 => 1064
# 16 => 2320
# 24 => 7664
# 32 => 17904
# 48 => 59728

# x = 2πa/λ
# |m|kd = 1
# k = 2π/λ

a = 5
λ = 2π

# m = 0.96+1.01im
# d = 1/(abs(m))

# m = 1.33 + .01im
# m = 2 + im

Nx = Ny = Nz = 8
origin = Point3(0., 0., 0.)
spacing =  [1., 1., 1.]


g = CartesianGrid(origin, spacing, (Nx, Ny, Nz))
center = DDA.center(g)
radius = (minimum(DDA.width(g)) + 1.49)/2


t = DDA.Sphere(center, radius)
r2, occ = DDA.dipoles(g, t)
length(r2)




Nx = Ny = Nz = 8
origin = Point3(-(Nx-1)/2, -(Ny-1)/2, -(Nz-1)/2)
spacing =  [1., 1., 1.]


g = CartesianGrid(origin, spacing, (Nx, Ny, Nz))
center = DDA.center(g)
radius = (minimum(DDA.width(g)) + 1.49)/2


t = DDA.Sphere(center, radius)
r2, occ = DDA.dipoles(g, t)
length(r2)




# abs(m)* k * d ≪ 1
# d * k  = 1 / (abs(m) )


############################################################################
# "Medium"" level interface to DDA calulations
############################################################################

# everything stored as a list and there is a mapping between idicies
# r (3, Nd)
# E (3, Nd)
# α (3, Nd)
# occ (3, Nd)
# A (Nx, Ny, Nz, 6) ??? lazy?

# using Rotations: RotZY
#
# R = RotZY(θ, ϕ)
# E₀ = R[:,1:2] * e
# kvec = R[:,3] * k


function calc(d)

    # parameters
    λ = 1.      # wavelength
    k = 2π / λ  # wavenumber
    # a = 1       # effective radii
    # d = .5      # spacing
    Nd = 8
    e = [1, 0]  # Jones polarisation vector
    θ, ϕ = 0, 0 # rotation angles [rad]

    # steps:
    # 1. create the coordinates of the dipoles,
    # CartesianGrid{Float64,3}(SVector(0.,0.,0.),SVector(1.,1.,1.),Dims((11,11,11)))
    # CartesianGrid{Float64,3}((0.,0.,0.),(1.,1.,1.),(11,11,11))

    Nx = Ny = Nz = Nd
    origin = Point3(0., 0., 0.)

    g = CartesianGrid(origin, [d, d, d], (Nx, Ny, Nz))

    center = DDA.center(g)
    radius = (minimum(DDA.width(g)) + 1.49*d)/2

    t = DDA.Sphere(center, radius)

    r2, occ = DDA.dipoles(g, t)
    N = length(r2)

    # a = effective radii ???
    a_eff = (3*N*d^3/(4*pi))^(1/3)


    # @btime DDA.dipoles($g, $t)


    # 2. assign the polarizability αj to each dipole

    m = 1.33 + .01im
    # m = 2 + im

    ε = m^2

    # ε = permitivity(Material.load("Ag.yaml"), k)

    # d = a * spacing # ???
    kvec = [0, 0, k]
    E₀ = [e..., 0]
    α = LDR(ε, d, kvec, E₀)

    alphas = fill(α, length(r2))


    # 3. calculated the incident field E_inc, at each dipole,
    pw = DDA.PlaneWave2(k, e, θ, ϕ)
    E_inc = field(pw, r2)


    # 4. assemble the interaction matrix A and
    A = DDA.interactions(k, r2, alphas)

    # 5. solve for P in the system of linear equations

    Ei = reinterpret(ComplexF64, E_inc)
    Pi = A \ Ei
    Pi = bicgstabl(A, Ei; abstol=1e-3, verbose=true)
    # Pi = cg(A, Ei; abstol=1e-3, maxiter=1000, verbose=true)
    # Pi = minres(A, Ei; abstol=1e-3, maxiter=1000, verbose=true)

    P = reinterpret(SVector{3,ComplexF64},Pi)

    # # The solution using FFT
    # A = InteractionTensor(g, k, α, occ)
    # P = bicgstabl(factorise(A), E; abstol=1e-8)


    # Analyse the results

    # @show DDA.C_abs(k, E₀, P, alphas)
    # @show DDA.C_ext(k, E₀, E_inc, P)
    # @show DDA.C_sca(k, E₀, E_inc, P, alphas)
    @show d
    return DDA.C_abs(k, E₀, P, alphas), DDA.C_ext(k, E₀, E_inc, P), DDA.C_sca(k, E₀, E_inc, P, alphas)

end

d = range(0, .5, length=101)[2:end]



R_abs = zeros(length(d))
R_ext = zeros(length(d))
R_sca = zeros(length(d))

for i in eachindex(d)
    r = calc(d[i])
    R_abs[i] = r[1]
    R_ext[i] = r[2]
    R_sca[i] = r[3]
end

a = d
Q_abs = @. R_abs / (π * a^2)
Q_ext = @. R_ext / (π * a^2)
Q_sca = @. R_sca / (π * a^2)


plt = plot(d, R_sca,
    # xlabel="ka",
    # ylabel="Scattering efficiency",
    # legend=false,
    # label="LABEL",
    # title="TITLE",
    yscale=:log10,
    ylims=(0.00005, maximum(R_sca))
)
plot!(d, R_abs)
display(plt)

plt = plot(d, R_sca./(π * d.^2),
# xlabel="ka",
# ylabel="Scattering efficiency",
# legend=false,
# label="LABEL",
# title="TITLE",
yscale=:log10,
ylims=(0.00005, maximum(R_sca./(π * d.^2))))
plot!(d, R_abs./(π * d.^2))
display(plt)


# calculate the field at a point
r = [1, 1, 1]
E = E_sca_FF(k, r, P, r_E)




############################################################################
# High-level interface for DDA calculations
############################################################################

# 1. Define a grid
Nx, Ny, Nz = 32, 32, 32
origin = [0, 0, 0]
spacing = [1, 1, 1]

g = CartesianGrid(origin, spacing, (Nx, Ny, Nz))

# 2. Define the target(s)
center = center(g)
radius = 5.

target = Sphere(center, radius)

# 3. Define the material properties
ε = 1.33 + 0.1im
α = LDR(ε)
scatterer = Scatterer(t, α)

# 4. Define incindent field
k = 2π      # wavenumber
e = [1, 0]  # Jones polarisation vector
θ, ϕ = 0, 0 # rotation angles [rad]

E_inc = PlaneWave(k, e, θ, ϕ)

# 5. Define the DDA problem
prob = DDAGridProblem(grid, scatterer, E_inc)

# 6. Solve the DDA problem
sol = solve(prob, BiCGStablFFT(), tol=1e-12)



# 3. Define the material properties
# TODO: wavelelngth or frequency???
# TODO: refractive index or permittivity???

Ag = get_material("main","Ag","Johnson")

λ_Ag = Ag.λ
ñ_Ag = Ag.n + Ag.k * im
ε_Ag = ñ_Ag.^2;
ω_Ag = ustrip.(2π * c_0 ./ (λ_Ag * 1u"μm") .|> u"THz")

ε = PermittivityTable(ω_Ag, ε_Ag)



# ε = Material.gold("") ???
# ε = Material.load("Ag.yaml")
#
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