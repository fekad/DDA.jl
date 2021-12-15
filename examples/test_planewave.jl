using Revise
using DDA


k = 1
e = [1, .3im]
θ = π/3
ϕ = π/5

pw1 = PlaneWave(k, e, θ, ϕ)
pw2 = DDA.PlaneWave2(k, e, θ, ϕ)

# r = rand(3, 10000)
#

r = [rand(3) for _ in 1:1000000]
f1(pw1::PlaneWave, r) = [field(pw1, r) for r in r]
f2(pw2::DDA.PlaneWave2, r) = [field(pw2, r) for r in r]

f1(pw1, r)
@time f1(pw1, r)
@time f2(pw2, r)