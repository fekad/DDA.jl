
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
    return 3d^3 / 4pi * (eps - 1) / (eps + 2)
end
# alpha_CM(eps, d) = 3d^3 / 4pi * (eps - 1) / (eps +2)


"""
LDR: Lattice Dispersion Relation

d: dipole "lattice" spacing
eps: complex permittivity
k_vec: wavevectors
E0: complex incident field amplitude

https://doi.org/10.1086/172396
https://doi.org/10.1364/JOSAA.11.001491
"""
function polarizability_LDR(eps, d, k, E₀)
    b1 = -1.8915316
    b2 =  0.1648469
    b3 = -1.7700004

    # wavenumber
    k_norm = norm(k)

    # a_hat: incident direction (unit vector)
    a_hat = normalize(k)

    # e_hat: polarisation state (unit vector)
    e_hat = normalize(E₀)

    S = sum((a_hat[i] * e_hat[i])^2 for i in 1:3)

    # First compute Clausius-Mossotti polarizability:
    alpha_CM = polarizability_CM(d, eps)

    # Determine polarizability by requiring that infinite lattice of
    # dipoles have dipersion relation of continuum.
    alpha_inf = alpha_CM / (1 + alpha_CM / d^3 * (b1 + eps * (b2 + b3 * S)) * (k_norm * d)^2)

    # Radiative-reaction correction:
    alpha = alpha_inf / (1 -  alpha_inf * 2 / 3 * im * k_norm^3 )

    return alpha
end



# """
# LDR: Lattice Dispersion Relation
#
# d: dipole "lattice" spacing
# eps: complex permittivity
# k_vec: wavevectors
# E0: complex incident field amplitude
# """
# function polarizability_LDR(d, eps, kvec, E0)
#     b1 = -1.8915316
#     b2 =  0.1648469
#     b3 = -1.7700004
#
#     # wavenumber
#     k = norm(kvec)
#
#     # a_hat: incident direction (unit vector)
#     a_hat = normalize(kvec)
#
#     # e_hat: polarisation state (unit vector)
#     e_hat = normalize(E0)
#
#     S = sum(@. (a_hat * e_hat)^2)
#
#     # First compute Clausius-Mossotti polarizability:
#     alpha_CM = polarizability_CM(d, eps)
#
#     return alpha_CM / (1 + alpha_CM / d^3 * ((b1 + eps * b2 + eps * b3 * S) * (k * d)^2 - 2 / 3 * im * (k * d)^3))
#
# end
#


