using DDA, Test

m = 1.33 + .1im  # imag. component to demonstrate absorption
eps = m^2
k = 2 * pi
d = 1 / (abs(m) * k)

alph = DDA.polarizability_CM(d, eps)
@test isapprox(alph, 8.351021041502063e-5 + 2.279588323948128e-5im; rtol=10^-4)

E0 = [1, 1, 0]
kvec = [0, 0, k]
alph = DDA.polarizability_LDR93(d, eps, kvec, E0)
@test isapprox(alph, 8.639466684174675e-5 + 2.5988440242162623e-5im; rtol=10^-4)

