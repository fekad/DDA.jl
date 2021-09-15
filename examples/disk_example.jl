using DDA

using plot

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


# incident field
E0 = 1
k = [0, 0, 1]
e = [1, 0, 0]
E_inc = E_inc()


# grid and target
r = 22
d = 2*r/30

s = Sphere(r, [0,0,0])
g = CubicGrid(-r:d:r,-r:d:r,-r:d:r)

# get_dipoles
# dipoles = get_dipoles(g, s)
dipoles = positions(g, s)

# material model

# solve
lambda = 550.
# lambda = 400:10:600
solve(dipoles, E_inc)

??material