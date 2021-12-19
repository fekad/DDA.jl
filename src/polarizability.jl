#
# !       CALPHA = polarizability prescription
# !              = 'LATTDR' for LDR of Draine & Goodman (1993)
# !              = 'GKDLDR' for LDR of Gutkowicz-Krusin & Draine (2004)
# !              = 'FLTRCD' for filtered coupled dipole approach
# !                             of Gay-Balmaz & Martin (2002) and
# !                             Yurkin, Min & Hoekstra (2010)
# !              [= 'SCLDR' not supported in present version]



# abstract type AbstractPolarizability end

# struct IsotropicPolarizability{T<: Number} <: Number
#     alpha::T
# end
#
# struct AnisotropicPolarizability{T<: Number} <: AbstractArray
#     alpha::T
# end

# isotropic (scalar or diagonal)
# anisotropic (symmetric tensor)
# homogenous vs. different

@doc raw"""
    CM(eps, d)

Calcualtes Clausius-Mossoti Polarizability of dipole array according to their refractive indexes `m` and lattice spacing `d`.

The Clausius-Mossoti polarizability of dipoles:
```math
\alpha_j^{CM} = \frac{3d^3}{4\pi} \frac{\varepsilon_j - 1}{\varepsilon_j + 2}
```

# Arguments
- `d`: Dipole lattice spacing
- `eps`: complex permitivirty of dipole
"""
CM(eps, d) = 3d^3 / 4π * (eps - 1) / (eps +2)

"""
    RR_correction(alpha, k)

Radiative-reaction correction
"""
RR_correction(alpha, k) = alpha / (1 -  alpha * 2 / 3 * im * k^3)


"""
    CMRR(eps, d, k)

"Clausius-Mossotti plus Radiative-Reaction" (CMRR)
See: B.T. Draine and J. Goodman, Astrophys. J. 405, pp. 685-697 (1993)
"""
function CMRR(eps, d, k)

    # First compute Clausius-Mossotti polarizability:
    alpha_CM = CM(eps, d)

    # Radiative-reaction correction:
    alpha = RR_correction(alpha_CM, k)

    return alpha

end


"""
    GOHG(eps, d, k)

Goedecke & O'Brien (1988) and Hage & Greenberg (1990) (GOHG)
See: B.T. Draine and J. Goodman, Astrophys. J. 405, pp. 685-697 (1993)
"""
function GOHG(eps, d, k)

    alpha_CM = CM(eps, d)

    alpha_NR =  alpha_CM / (1 - (4π/3)^(1/3) * alpha_CM * k^2 / d)

    # Radiative-reaction correction:
    alpha = RR_correction(alpha_NR, k)

    return alpha

end


"""
    ILDR(eps, d, k)

"Isotropized" Lattice Dispersion Relation (ILDR)
See: B.T. Draine and J. Goodman, Astrophys. J. 405, pp. 685-697 (1993)
"""
function ILDR(eps, d, k)
    b1 = -1.8915316
    b2 = -0.1891532

    # First compute Clausius-Mossotti polarizability:
    alpha_CM = CM(eps, d)

    alpha_NR = alpha_CM / (1 + alpha_CM / (b1 + eps * b2) * k^2 / d)

    # Radiative-reaction correction:
    alpha = RR_correction(alpha_NR, k)

    return alpha

end


"""
    LDR(eps, d, kvec, E₀)

Lattice Dispersion Relation (LDR)

- https://doi.org/10.1086/172396
- https://doi.org/10.1364/JOSAA.11.001491

# Arguments
- `d`: dipole "lattice" spacing
- `eps`: complex permittivity
- `k_vec`: wavevectors
- `E0`: complex incident field amplitude

"""
function LDR(eps, d, kvec, E₀)
    b1 = -1.8915316
    b2 =  0.1648469
    b3 = -1.7700004

    # wavenumber
    k = norm(kvec)

    # a_hat: incident direction (unit vector)
    a_hat = normalize(kvec)

    # e_hat: polarisation state (unit vector)
    e_hat = normalize(E₀)

    S = sum(1:3) do i
        (a_hat[i] * e_hat[i])^2
    end
    # S = sum((a_hat[i] * e_hat[i])^2 for i in 1:3)

    # First compute Clausius-Mossotti polarizability:
    alpha_CM = CM(eps, d)

    # Determine polarizability by requiring that infinite lattice of
    # dipoles have dipersion relation of continuum.
    alpha_inf = alpha_CM / (1 + alpha_CM * (b1 + eps * (b2 + b3 * S)) * k^2 / d)

    # Radiative-reaction correction:
    alpha = RR_correction(alpha_inf, k)

    @assert alpha ≈ alpha_CM / (1 + alpha_CM / d^3 * ((b1 + eps * b2 + eps * b3 * S) * (k * d)^2 - 2 / 3 * im * (k * d)^3))

    return alpha

end

