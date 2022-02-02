
struct Atomic <: AbstractPolarizability
    ε
end


struct LDRModel <: AbstractPolarizability
    ε
end


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

Calculates Clausius-Mossoti polarizability for pointlike atoms (dipoles) on a cubic lattice where the lattice spacing is `d`.

The Clausius-Mossoti polarizability of dipoles:
```math
\alpha_j^{CM} = d^3 \frac{3}{4\pi} \frac{\varepsilon_j - 1}{\varepsilon_j + 2}
```

# Arguments
- `eps`: complex permitivirty of dipole
- `d`: Dipole lattice spacing
"""
CM(eps, d) = d^3 * 3 / 4π * (eps - 1) / (eps + 2)

"""
    RR_correction(alpha, k)

"Radiative-reaction" correction (according to the optical theorem)
See: B.T. Draine and J. Goodman, Astrophys. J. 405, pp. 685-697 (1993)
"""
# RR_correction(alpha, k, d) = alpha / d^3 / (1 - alpha * 2 / 3 * im * (k * d)^3)
RR_correction(alpha, k) = alpha / (1 - alpha * 2 / 3 * im * k^3)



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

    alpha_NR = alpha_CM / (1 - (4π / 3)^(1 / 3) * alpha_CM * k^2 / d)

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
# function LDR(eps, d, k, knorm, Enorm)
function LDR(eps, d, kvec, E₀)
    b1 = -1.8915316
    b2 = 0.1648469
    b3 = -1.7700004

    # wavenumber
    k = norm(kvec)

    # a_hat: incident direction (unit vector)
    a_hat = normalize(kvec)

    # e_hat: polarisation state (unit vector)
    e_hat = normalize(E₀)

    S = sum(1:3) do i
        # (a_hat[i] * abs(e_hat[i]))^2
        (a_hat[i] * e_hat[i])^2
    end
    # S = sum((a_hat[i] * e_hat[i])^2 for i in 1:3)

    # First compute Clausius-Mossotti polarizability:
    alpha_CM = CM(eps, d)

    # Determine polarizability by requiring that infinite lattice of
    # dipoles have dipersion relation of continuum.
    alpha_inf = alpha_CM / (1 + alpha_CM / d^3 * (b1 + eps * (b2 + b3 * S)) * (k * d)^2)

    # Radiative-reaction correction:
    alpha = RR_correction(alpha_inf, k)

    @assert alpha ≈ alpha_CM / (1 + alpha_CM / d^3 * ((b1 + eps * b2 + eps * b3 * S) * (k * d)^2 - 2 / 3 * im * (k * d)^3))

    return alpha

end


#
# GKDLDR
# ! B1 = (c_1/pi)*ak2 = (-5.9424219/pi)*ak2 = -1.8915316*ak2
# ! B2 = (c_2/pi)*ak2 = (0.5178819/pi)*ak2 = 0.1648469*ak2
# ! B3 = -[(3c_2+c_3)/pi]*ak2 = -[(3*0.5178819+4.0069747)/pi]*ak2
# !                           = -1.7700004*ak2
# ! B3L = B3*A(I)**2 where a_i = unit vector in direction of propagation
#
#          B1=-1.8915316_WP*AK2
#          B2=0.1648469_WP*AK2
#          B3=-1.7700004_WP*AK2
#          DO L=1,3
#             B3L=B3*AK_TF(L)*AK_TF(L)/AK2
#             DO IA=1,NAT
#                IC=ICOMP(IA,L)
#                IF(IC>0)THEN
#
# !*** First compute Clausius-Mossotti polarizability:
#
#                   CXTERM=(.75_WP/PI)*(CXEPS(IC)-1._WP)/(CXEPS(IC)+2._WP)
#
# !*** Determine polarizability by requiring that infinite lattice of
# !    dipoles have dipersion relation of continuum.
#
#
#                   CXTERM=CXTERM/(1._WP+CXTERM*(B1+CXEPS(IC)*(B2+B3L)))
#
# !*** Radiative-reaction correction:
#
#                   CXALPH(IA,L)=CXTERM/(1._WP+CXTERM*CXRR)
#
# ! set off-diagonal terms to zero
#
#                   CXALOF(IA,L)=0._WP
#
#                ELSEIF(IC==0)THEN
#
# ! To avoid divisions by zero, etc., set CXALPH=1 for vacuum sites.
#
#                   CXALPH(IA,L)=1._WP
#                   CXALOF(IA,L)=0._WP
#                ENDIF
#             ENDDO
#          ENDDO

# function FCD_Pol(DF,dV,dx,i_energy)
#     !
#     ! "Filtered Coupled Dipole Method" for very large refractive indices
#     ! See: Yurkin, M. A., Min, M., & Hoekstra, A. G. (2010).
#     ! Application of the discrete dipole approximation to very large refractive indices:
#     ! Filtered coupled dipoles revived. Physical Review E, 82(3), 036703
#     !
#     use kinds             , only: pr
#     use physical_constants, only: pi, eV_to_k0
#     use inputs            , only: EnLoss
#
#     implicit none
#     integer(kind=pr)                , intent(in) :: i_energy
#     complex(kind=pr)                , intent(in) :: DF
#     complex(kind=pr)                             :: CM, D, FCD_Pol, cst1, cst2, cst3
#     real   (kind=pr)              , dimension(3) :: k
#     real   (kind=pr)                , intent(in) :: dV, dx
#
#     CM = CM_Pol(DF,dV)
#
#     ! wave number
#     k(1) = EnLoss(i_energy)*eV_to_k0
#     k(2) = k(1)*k(1) ; k(3) = k(2)*k(1)
#
#     cst1 = (1.3333333333333333_pr)*k(2)*(dx*dx)
#
#     ! 4/3 = 1.3333333333333333
#
#     cst2 = (0.21220658488400759_pr)*log((pi - k(1)*dx)/(pi + k(1)*dx))*k(3)*dV
#
#     ! 2/(3*pi) = 0.21220658488400759
#
#     cst3 = (0.66666666666666663_pr)*(cmplx(0.0_pr,1.0_pr))*k(3)*dV
#
#     ! 2/3 = 0.66666666666666663
#
#     D = (CM/dV)*(cst1 + cst2 + cst3)
#
#     FCD_Pol = CM/(1.0_pr + D)
#
#     return
#
# end function FCD_Pol