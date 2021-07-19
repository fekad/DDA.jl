using Revise
using DDA, Test


k = 6.283185307179586
E0 =[1, 1, 0]
kvec = [0, 0, k]

r = [
    [-0.0597, -0.0597, -0.0597],
    [-0.0597, -0.0597, 0.0597],
    [-0.0597, 0.0597, -0.0597],
    [-0.0597, 0.0597, 0.0597],
    [0.0597, -0.0597, -0.0597],
    [0.0597, -0.0597, 0.0597],
    [0.0597, 0.0597, -0.0597],
    [0.0597, 0.0597, 0.0597]
]

Ei_exp = [
    0.9306 - 0.3662im,
    0.9306 - 0.3662im,
    0.0000 + 0.0000im,
    0.9306 + 0.3662im,
    0.9306 + 0.3662im,
    0.0000 + 0.0000im,
    0.9306 - 0.3662im,
    0.9306 - 0.3662im,
    0.0000 + 0.0000im,
    0.9306 + 0.3662im,
    0.9306 + 0.3662im,
    0.0000 + 0.0000im,
    0.9306 - 0.3662im,
    0.9306 - 0.3662im,
    0.0000 + 0.0000im,
    0.9306 + 0.3662im,
    0.9306 + 0.3662im,
    0.0000 + 0.0000im,
    0.9306 - 0.3662im,
    0.9306 - 0.3662im,
    0.0000 + 0.0000im,
    0.9306 + 0.3662im,
    0.9306 + 0.3662im,
    0.0000 + 0.0000im
]


Ei = E_inc(E0, kvec, r)
@test isapprox(Ei, Ei_exp; rtol=1e-3)