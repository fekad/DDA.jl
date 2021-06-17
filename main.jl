using Revise
using DDA
using Plots

k = [0, 0, 1.]


E_0 = 1.

E_inc = DDA.PlaneWave(E_0, k)


r_n = [1, 1, 1]
r = r_n * (0:0.1:10)'

E = DDA.field(r[:,1], E_inc)
E = DDA.field(r, E_inc)


plot(real(E))


