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


g = CartesianGrid([.5,.5,.5], ones(3), (30, 30, 30))
s = Sphere(15., [15, 15, 15])
d = dipoles(g, s)

dd = reduce(hcat, d)
scatter3d(dd[1,:],dd[2,:],dd[3,:])

# g = CartesianGrid([-14.5, -14.5, -14.5], ones(3), (30,30,30))
# s = Sphere(15., [0, 0, 0])
# d = dipoles(g,s)

g = CartesianGrid([.5,.5,.5], ones(3), (30, 30, 30))
s = Disk(15, 5, [15, 15, 0 ])
d = dipoles(g, s)
dd = reduce(hcat, d)

# gr()
plotlyjs()
scatter(dd[1,:],dd[2,:],dd[3,:], aspect_ratio=:auto)
xlims!(0,30); ylims!(0, 30); zlims!(0, 30)

