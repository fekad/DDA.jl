
################################################################
# 1. Define incident field

k = 2π          # wave number
E0 = [1., 1, 0]
kvec = [0, 0, k]

E_inc(E0, kvec, e)

################################################################
# 2. Define grid

g = CubicGrid(origin, spacing, dims)

################################################################
# 3. Define structure (Define material model?)
# s = Sphere(R, [0, 0, 0])
# S = dipoles(g, s)
# alph = polarizability_LDR(m, d, kvec, E0)
# OR
# s = Sphere(R, [0, 0, 0], eps)
# S = dipoles(s, g, polarisability_model=:LDR) ??? d, kvec, E0
# OR
# s = Sphere(R, [0, 0, 0], eps, polarisability_model=:LDR) ??? d, kvec, E0
# S = dipoles(s, g)
#
# Composition?????

s = Sphere(R, [0, 0, 0])

S = dipoles(g, s)
N = length(S)

ε = m^2
polarizability_CM(d, eps)
polarizability_LDR(mi, d, kvec, E0)
alph = polarizability_LDR(eps, d, kvec, E₀)

# list of dipoles + polarisability OR grid + shape + polarisability model
prob = DDAProblem(s, g, polarisability_model=:LDR, E_inc)


################################################################
# 4. Solve the linear problem
# Simple()
# Gmres()
# Minres()
# Qmr()
# struct SimpleTsit5 <: DiffEqBase.AbstractODEAlgorithm end
# export SimpleTsit5

# Solver
# 1. create the coordinates of the dipoles,
# 2. assign the polarizability αj to each dipole,
# 3. calculated the incident field E_inc, at each dipole,
# 4. assemble the interaction matrix A and
# 5. solve for P in the system of linear equations
sol = solve(prob, Simple(), reltol=1e-8, abstol=1e-8, maxiter)


C_ext(sol)
C_abs(sol)
C_sca(sol)








# Low-level inteface

# 1. create the coordinates of the dipoles,
g = CubicGrid(origin, spacing, dims)
s = Sphere(radius, center)
dipoles(g, s)

# 2. assign the polarizability αj to each dipole,

# 3. calculated the incident field E_inc, at each dipole,
Ei = E_inc(E0, kvec, r)

# 4. assemble the interaction matrix A and
A = interaction(k, dipoles, alphas)
calc_Ajj(α)
calc_Ajk(k, rj, rk)

# 5. solve for P in the system of linear equations
solve(A, E_incs)
# P = A \ Ei
# P = gmres(A, Ei)
# P = minres(A, Ei)
# P = qmr(A, Ei)


C_abs(k, E0, P, alph)
C_ext(k, E0, Ei, P)
C_sca(k, E0, Ei, P, alph)

E_sca_FF(k, r, P, r_E)




