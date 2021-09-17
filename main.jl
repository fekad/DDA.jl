using Revise
using DDA
using Plots

plotlyjs()

# k = [0, 0, 1.]
#
#
# E_0 = 1.
#
# E_inc = DDA.PlaneWave(E_0, k)
#
#
# r_n = [1, 1, 1]
# r = r_n * (0:0.1:10)'
#
# E = DDA.field(r[:,1], E_inc)
# E = DDA.field(r, E_inc)
#
#
# plot(real(E))


# target: Au disk
# d = 10 nm
# # of dipoles along the diameter: 50
# grid spacing: d = 0.2 nm
# # of dipoles: ??
# spaceing: d = (1) 0.2 nm, (2) 0.4 nm, (5) 1 nm,  (100) 20 nm
# >> peak resonance at ??? nm

# steps:
# 1. create the coordinates of the dipoles,
# 2. assign the polarizability αj to each dipole,
# 3. calculated the incident field E_inc, at each dipole,
# 4. assemble the interaction matrix A and
# 5. solve for P in the system of linear equations




# 1. create the coordinates of the dipoles,

g = CartesianGrid([1.5,1.5,1.5], 3*ones(3), (30, 30, 30))
s = Sphere(15., [15, 15, 15])
r = dipoles(g, s)

rr = reduce(hcat, r)
scatter3d(rr[1,:],rr[2,:],rr[3,:])

# g = CartesianGrid([-14.5, -14.5, -14.5], ones(3), (30,30,30))
# s = Sphere(15., [0, 0, 0])
# d = dipoles(g,s)



g = CartesianGrid([-100, -100, -100], 200*ones(3), (11, 11, 11))
s = Disk(15, 5, [15, 15, 0 ])
r = dipoles(g, s)
rr = reduce(hcat, r)



g = CartesianGrid([-200, -200, -200], 40*ones(3), (11, 11, 11))
s = Sphere(200., [0, 0, 0])
r = dipoles(g, s)
rr = reduce(hcat, r)


# gr()
plotlyjs()
scatter(rr[1,:],rr[2,:],rr[3,:], aspect_ratio=:auto)
xlims!(-200, 200); ylims!(-200, 200); zlims!(-200, 200)



# 2. assign the polarizability αj to each dipole,

E₀ = [1, 0, 0]
k = 2π * [0, 0, 1]

n = 4.686878048780487+0.08797560975609756im
eps = n^2
d = g.spacing[1]


alpha = polarizability_LDR(eps, d, k, E₀)


# 3. calculated the incident field E_inc, at each dipole,

p = PlaneWave(E₀, k)
E_inc = field(r, p)


# 4. assemble the interaction matrix A and
using LinearAlgebra

@time A = DDA.interaction(norm(k), r, alpha * ones(length(r)))



# 5. solve for P in the system of linear equations


using IterativeSolvers

@time s1=cg(A, E_inc[:]; abstol=1e-8, maxiter=1000, verbose=true)
@time s2=minres(A, E_inc[:]; abstol=1e-8, maxiter=10000)
@time s3=bicgstabl(A, E_inc[:]; abstol=1e-8)

@time s4=A \ E_inc[:]