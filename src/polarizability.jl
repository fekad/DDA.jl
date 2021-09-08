
"""
polarizability_CM(eps, d)

Calcualtes Clausius-Mossoti Polarizability of dipole array according to their refractive indexes `m` and lattice spacing `d`.

The Clausius-Mossoti polarizability of dipoles:
```math
\\alpha_j^{CM} = \\frac{3d^3}{4\\pi} \\frac{\\varepsilon_j - 1}{\\varepsilon_j + 2}
```

# Arguments:
- `d`: Dipole lattice spacing
- `eps`: complex permitivirty of dipole
"""
function polarizability_CM(d, eps)
    return 3d^3 / 4pi * (eps - 1) / (eps +2)
end
# alpha_CM(eps, d) = 3d^3 / 4pi * (eps - 1) / (eps +2)



function polarizability_LDR93(d, eps, kvec, E0)
    b1 = -1.8915316
    b2 =  0.1648469
    b3 = -1.7700004

    k = norm(kvec)
    a_hat = normalize(kvec)
    e_hat = normalize(E0)
    S = sum(@. (a_hat * e_hat)^2)

    # First compute Clausius-Mossotti polarizability:
    alpha_CM = polarizability_CM(d, eps)

    return alpha_CM / (1 + alpha_CM / d^3 * ((b1 + eps * b2 + eps * b3 * S) * (k * d)^2 - 2/3 * im * (k * d)^3))

end



# function polarizability_LDR93_5(eps, d, kvec, E0)
#     b1 = -1.8915316
#     b2 =  0.1648469
#     b3 = -1.7700004
#
#     k = norm(kvec)
#     a_hat = normalize(kvec)
#     e_hat = normalize(E0)
#     S = sum(@. (a_hat * e_hat)^2)
#
#     # First compute Clausius-Mossotti polarizability:
#     alpha_CM = polarizability_CM(d, eps)
#
#     # Determine polarizability by requiring that infinite lattice of
#     # dipoles have dipersion relation of continuum.
#     alpha_inf = alpha_CM / (1 + alpha_CM / d^3 * (b1 + eps * (b2 + b3 * S)) * (k * d)^2)
#
#     # Radiative-reaction correction:
#     alpha = k^3 * alpha_inf/(1 -  2/3 * im * k^3* alpha_inf)
#     #   NOTE: k^3 ???? d^3 ????
#     return alpha
#
# end
# alph = polarizability_LDR93_5(d, eps, kvec, E0)
#
# function polarizability_LDR93_4(d, eps, kvec, E0)
#     b1 = -1.8915316
#     b2 =  0.1648469
#     b3 = -1.7700004
#
#     k = norm(kvec)
#     a_hat = normalize(kvec)
#     e_hat = normalize(E0)
#     S = sum(@. (a_hat * e_hat)^2)
#
#     # First compute Clausius-Mossotti polarizability:
#     alpha_CM = polarizability_CM(d, eps)
#
#     # Determine polarizability by requiring that infinite lattice of
#     # dipoles have dipersion relation of continuum.
#     alpha_inf = alpha_CM / (1 + alpha_CM/ d^3  * (b1 + eps * ( b2 + b3 * S)) * (k * d)^2)
#
#     # Radiative-reaction correction:
#     alpha =  alpha_inf/(1 - 2/3 * im * (k * d)^3 * alpha_inf)
#
#     return alpha
#
#
# end
# alph = polarizability_LDR93_4(d, eps, kvec, E0)





