using DDA, Test


k = 6.283185307179586
e = [1, 1]

pw = PlaneWave(k, e, 0, 0)

r = [-0.0597, -0.0597, -0.0597]
Ei = field(pw, r)
Ei_exp = [
    0.9306 - 0.3662im,
    0.9306 - 0.3662im,
    0.0000 + 0.0000im,
]
@test isapprox(Ei, Ei_exp; rtol = 1e-3)


r = [-0.0597, -0.0597, 0.0597]
Ei= field(pw, r)
Ei_exp = [
    0.9306 + 0.3662im,
    0.9306 + 0.3662im,
    0.0000 + 0.0000im
]
@test isapprox(Ei, Ei_exp; rtol = 1e-3)

r = [-0.0597, 0.0597, -0.0597]
Ei= field(pw, r)
Ei_exp = [
    0.9306 - 0.3662im,
    0.9306 - 0.3662im,
    0.0000 + 0.0000im
]
@test isapprox(Ei, Ei_exp; rtol = 1e-3)

r = [-0.0597, 0.0597, 0.0597]
Ei= field(pw, r)
Ei_exp = [
    0.9306 + 0.3662im,
    0.9306 + 0.3662im,
    0.0000 + 0.0000im
]
@test isapprox(Ei, Ei_exp; rtol = 1e-3)

r = [0.0597, -0.0597, -0.0597]
Ei= field(pw, r)
Ei_exp = [
    0.9306 - 0.3662im,
    0.9306 - 0.3662im,
    0.0000 + 0.0000im
]
@test isapprox(Ei, Ei_exp; rtol = 1e-3)

r = [0.0597, -0.0597, 0.0597]
Ei= field(pw, r)
Ei_exp = [
    0.9306 + 0.3662im,
    0.9306 + 0.3662im,
    0.0000 + 0.0000im
]
@test isapprox(Ei, Ei_exp; rtol = 1e-3)

r = [0.0597, 0.0597, -0.0597]
Ei= field(pw, r)
Ei_exp = [
    0.9306 - 0.3662im,
    0.9306 - 0.3662im,
    0.0000 + 0.0000im
]
@test isapprox(Ei, Ei_exp; rtol = 1e-3)

r = [0.0597, 0.0597, 0.0597]
Ei= field(pw, r)
Ei_exp = [
    0.9306 + 0.3662im,
    0.9306 + 0.3662im,
    0.0000 + 0.0000im
]
@test isapprox(Ei, Ei_exp; rtol = 1e-3)
