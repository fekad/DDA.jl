using Revise
using DDA

using Plots
plotlyjs()


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
d = 2r / 15

s = Sphere(r, [0,0,0])
# @. ((1:30) - .5) * d - r ≈ -r + d/2:d:r - d/2
g = CubicGrid(-r + d / 2:d:r - d / 2, -r + d / 2:d:r - d / 2, -r + d / 2:d:r - d / 2)


# get_dipoles
# dipoles = get_dipoles(g, s)
dipoles = positions(g, s)

# material model

# solve
lambda = 550.
# lambda = 400:10:600
solve(dipoles, E_inc)









E0 = [1., 1, 0]

m1 = 1.33 # relative refractive index of water
# m1 = 1.33 + .1im   # imag. component to demonstrate absorption
k = 2π          # wave number
d = 1 / (abs(m1) * k) # lattice spacing


# number of dipoles in the approximate sphere; more is required as the
# radius increases
# nrange = [8 32 136 280 552 912 1472]
# nrange = [8 32 136 280 552 912 1472 2176 3112 4224 5616 7208 9328 11536]
g = CartesianGrid([-20.5,-20.5,-20.5], ones(3), (41, 41, 41))



aeff = zeros(length(nrange))
Cext = zeros(length(nrange))
Cabs = zeros(length(nrange))
Cscat = zeros(length(nrange))

Rrange = 1.:14.
# for (ix, R) = enumerate(Rrange)
(ix, R) = first(enumerate(Rrange))

# 1. calculated the incident field E_inc, at each dipole,

kvec = [0, 0, k]
Ei = E_inc(E0, kvec, r)

# 2. create the coordinates of the dipoles,

s = Sphere(R, [0, 0, 0])
S = dipoles(g, s)
N = length(S)

r = d * S # scaling

# 3. assign the polarizability αj to each dipole,

m = m1 * ones(N)
alph = [DDA.polarizability_LDR(mi, d, kvec, E0) for mi in m]

# 4. assemble the interaction matrix A and
using LinearAlgebra
A = DDA.interaction(norm(k), r, alph)

# 5. solve for P in the system of linear equations

P = A \ Ei
# P = gmres(A, Ei)
# P = minres(A,Ei)
# P = qmr(A,Ei)

# 5. collect results

# # the corresponding effective radii of the spheres
aeff[ix] = d * (3 * N / 4π).^(1 / 3)

Cext[ix] = DDA.C_ext(k, E0, Ei, P)
Cabs[ix] = DDA.C_abs(k, E0, Ei, P, alph)
Cscat[ix] = Cext[ix] - Cabs[ix]
# end


# # Here, we plot the efficiencies Q instead of the cross sections C
# # Q = C/(pi*r^2)
# plot(k*arange,Cext./(pi*arange.^2),'*')
# plot(k*arange,Cabs./(pi*arange.^2),'x')
# plot(k*arange,Cscat./(pi*arange.^2),'o')
# legend('Q_{ext}','Q_{abs}','Q_{scat}',0,0)
# ylabel('Q')
# xlabel('2\pia/\lambda') # size parameter ka
# title(['m = ' num2str(m1)])





nothing
# =============================================================================

#  phase function demo for a water sphere

# Parameters
k = 2π #  wave number
E0 = [1, 0, 0] #  x-polarization
m1 = 1.33
d = 1 / (abs(m1) * k)
a_eff = (3N / 4 * π)^(1 / 3) * 1 / (k * abs(m1))


# shapepath = "../../shape/"
# rfile = shapepath * "sphere_552.txt"
# S = dlmread(rfile)

g = CartesianGrid([-20.5,-20.5,-20.5], ones(3), (41, 41, 41))
s = Sphere(5., [0, 0, 0])
S = dipoles(g, s)

N = length(S)

r = d * S
m = m1 * ones(N)

kvec = [0, 0, k] # propagating in +z direction
# Ei = E_inc(E0, kvec, r)
E_inc = PlaneWave(E0, kvec)
Ei = field(r, E_inc)

alph = polarizability_LDR(d, m, k)
interaction_A(k,r,alph)
P = gmres(A, Ei)

theta_range = range(0, pi, length=100)

Esca_S = zeros(length(theta_range))
Esca_P = zeros(length(theta_range))
Einc_S = zeros(length(theta_range))
Einc_P = zeros(length(theta_range))
#
# for (ix, theta) = enumerate(theta_range)
#
#     phi = 90 #  perpendicular to x-z plane
#     r_E = zeros(1, 3) #  evaluation point
#     r_E[1], r_E[2], r_E[3] = rtp2xyz(100, theta, phi)
#     E = E_sca_FF(k, r, P, r_E)
#     Esca_S[ix] = norm(E)
#     kr = dot([k, k, k], r_E, 2)
#     expikr = exp(i .* kr)
#     E1 = [E0[1] * expikr ,E0[2] * expikr, E0[3] * expikr]
#     Einc_S[ix] = norm(E1)
#
#
#     phi = 0 #  parallel to x-z plane
#     r_E = zeros(1, 3)
#     [r_E(1) r_E(2) r_E(3)] = rtp2xyz(100, theta, phi)
#     E = E_sca_FF(k, r, P, r_E)
#     Esca_P[ix] = norm(E)
#     kr = dot([k, k, k], r_E, 2)
#     expikr = exp(i .* kr)
#     E1 = [E0(1) * expikr E0(2) * expikr E0(3) * expikr]
#     Einc_P[ix] = norm(E1)
# end
#
# ph = semilogy(range * 180 / pi, Esca_P.^2 ./ Einc_P.^2, "--", range * 180 / pi, Esca_S.^2 ./ Einc_S.^2)
# xlim([0 180])
# ylabel("log |E_{sca}|^2")
# xlabel("phase angle")
# title(["ka = " num2str(k * a_eff) ", m = " num2str(m1) ", N = " int2str(N)],"FontSize",14)
# legend("parallel","perpendicular")
#
#
