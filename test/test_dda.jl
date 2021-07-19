using Revise
using DDA, Test

# using LinearAlgebra
#
# function calc_Aj(k, r, alph, im, blockdiag)
# #     r_jk = norm(r)
# #     r_n = r / r_jk
# #
# #     A_jk = exp(im * k * r_jk) / r_jk * (k^2 * (r_n * r_n' - I) + (im * k * r_jk - 1) / (r_jk^2) * (3 * r_n * r_n' - I))
#
#     A_jk = zeros(ComplexF64, length(r), 3, 3)
#     for (k) in
#     r_jk = norm(r)
#     r_n = r / r_jk
#
#     A_jk = exp(im * k * r_jk) / r_jk * (k^2 * (r_n * r_n' - I) + (im * k * r_jk - 1) / (r_jk^2) * (3 * r_n * r_n' - I))
#     return A_jk
# end




alph = 1.0e-04 * [
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im]
k = 6.283185307179586
j = 0
blockdiag = 1

r = [
    [0, 0, 0],
    [1, 1, 1]
]

Aj_exp = [
    [10614.7014265967 - 3193.37990597579im, 0.00000000000000 + 0.00000000000000im,
        0.00000000000000 + 0.00000000000000im, 1.71006105686162 + 15.0987189646511im,
        1.24771766013322 - 7.59383349109870im, 1.24771766013322 - 7.59383349109870im],
    [0.00000000000000 + 0.00000000000000im, 10614.7014265967 - 3193.37990597579im,
        0.00000000000000 + 0.00000000000000im, 1.24771766013322 - 7.59383349109870im,
        1.71006105686162 + 15.0987189646511im, 1.24771766013322 - 7.59383349109870im],
    [0.00000000000000 + 0.00000000000000im, 0.00000000000000 + 0.00000000000000im,
        10614.7014265967 - 3193.37990597579im, 1.24771766013322 - 7.59383349109870im,
        1.24771766013322 - 7.59383349109870im, 1.71006105686162 + 15.0987189646511im],
]

Aj = calc_Aj(k, r, alph, im, blockdiag)

@test isapprox(Aim, Aj_exp; rtol=1e-3)



alph = 1.0e-04 * [
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im]
k = 6.283185307179586
j = 0
blockdiag = 1

r = [[0, 0, 0],
                    [1, 1, 1],
                    [-1, -1, -1]]

Aj_exp = [
    [10614.7014265967 - 3193.37990597579im, 0.00000000000000 + 0.00000000000000im,
        0.00000000000000 + 0.00000000000000im, 1.71006105686162 + 15.0987189646511im,
        1.24771766013322 - 7.59383349109870im, 1.24771766013322 - 7.59383349109870im,
        1.71006105686162 + 15.0987189646511im, 1.24771766013322 - 7.59383349109870im,
        1.24771766013322 - 7.59383349109870im],
    [0.00000000000000 + 0.00000000000000im, 10614.7014265967 - 3193.37990597579im,
        0.00000000000000 + 0.00000000000000im, 1.24771766013322 - 7.59383349109870im,
        1.71006105686162 + 15.0987189646511im, 1.24771766013322 - 7.59383349109870im,
        1.24771766013322 - 7.59383349109870im, 1.71006105686162 + 15.0987189646511im,
        1.24771766013322 - 7.59383349109870im],
    [0.00000000000000 + 0.00000000000000im, 0.00000000000000 + 0.00000000000000im,
        10614.7014265967 - 3193.37990597579im, 1.24771766013322 - 7.59383349109870im,
        1.24771766013322 - 7.59383349109870im, 1.71006105686162 + 15.0987189646511im,
        1.24771766013322 - 7.59383349109870im, 1.24771766013322 - 7.59383349109870im,
        1.71006105686162 + 15.0987189646511im],
]

Aj = calc_Aj(k, r, alph, im, blockdiag)


@test isapprox(Aim, Aj_exp; rtol=1e-3)






alph = 1.0e-04 * [
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im,
    0.8639 + 0.2599im]

k = 6.283185307179586
blockdiag = 1

r = [[0, 0, 0], [1, 1, 1]]

A_exp = [
    [10614.70 - 3193.380im, 0.000000 + 0.000000im, 0.000000 + 0.000000im, 1.710061 + 15.09872im,
        1.247718 - 7.593833im, 1.247718 - 7.593833im],
    [0.000000 + 0.000000im, 10614.70 - 3193.380im, 0.000000 + 0.000000im, 1.247718 - 7.593833im,
        1.710061 + 15.09872im, 1.247718 - 7.593833im],
    [0.000000 + 0.000000im, 0.000000 + 0.000000im, 10614.70 - 3193.380im, 1.247718 - 7.593833im,
        1.247718 - 7.593833im, 1.710061 + 15.09872im],
    [1.710061 + 15.09872im, 1.247718 - 7.593833im, 1.247718 - 7.593833im, 10614.70 - 3193.380im,
        0.000000 + 0.000000im, 0.000000 + 0.000000im],
    [1.247718 - 7.593833im, 1.710061 + 15.09872im, 1.247718 - 7.593833im, 0.000000 + 0.000000im,
        10614.70 - 3193.380im, 0.000000 + 0.000000im],
    [1.247718 - 7.593833im, 1.247718 - 7.593833im, 1.710061 + 15.09872im, 0.000000 + 0.000000im,
        0.000000 + 0.000000im, 10614.70 - 3193.380im]
]

A = interaction_A(k, r, alph, blockdiag)

@test isapprox(A, A_exp; rtol=1e-4)
