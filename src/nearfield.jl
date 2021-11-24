

@doc raw"""
In the far field, the scattered electric field is given by:
$$
E_{sca} = \frac{k^2 \exp(ikr)}{r} \sum \limits_{j=1}^{N} \exp(-ik\hat{r} \cdot r_j)(\hat{r}\hat{r}-1_3)P_j
$$
"""
function E_sca_FF(k, r, P, r_E)
    # k: wave number
    # r: dipole coordinates (N x 3 matrix)
    # P: polarizations (vector of length 3N; Px1,Py1,Pz1 ... PxN,PyN,PzN)
    # r_E: coord for the point at which to calculate the far field
    # Note: coordinates are relative to origin
    N = length(r)

    E = zeros(ComplexF64, 3)
    r_norm = norm(r_E)
    r_hat = r_E / r_norm
    A = k^2 * exp(im * k * r_norm) / r_norm
    for i in 1:N
        E += exp(-im * k * dot(r_hat, r[i])) * (r_hat' .* r_hat - I) * P[:,i]
    end
    return A * E
end