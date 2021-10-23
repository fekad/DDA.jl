using DDA

using Plots

plotlyjs()

# Target:
# r = 22 nm
# # of dipoles along the diameter: 30
# grid spacing: d = 44/30 nm
# # of dipoles: 14328
# >> peak resonance at 500 nm

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
d = 2r/15

s = Sphere(r, [0,0,0])
# @. ((1:30) - .5) * d - r ≈ -r + d/2:d:r - d/2
g = CubicGrid(-r+d/2:d:r-d/2,-r+d/2:d:r-d/2,-r+d/2:d:r-d/2)


# get_dipoles
# dipoles = get_dipoles(g, s)
dipoles = positions(g, s)

# material model

# solve
lambda = 550.
# lambda = 400:10:600
solve(dipoles, E_inc)

# ??material



