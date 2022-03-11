"""Notes:
Calculating far-field electric field:
- far-field on a shpere on the surface of a box???
- far-field fieture: radiation pattern?
- far-field electric field should have a origin...
- expand it into Spherical polynomials harmonics

Calculating near-field electric field:
- in genereal near field electric field only make sense outside of the target
- inside the target can only be calculated for grid

"""

@doc raw"""
    E_sca(k, coords, P, r)

In the near field, the scattered electric field is given by:
```math
E_{sca} = \sum \limits_{j=1}^{N} A(r_j,r_E) P_j
```

Note: it is divergent at the position of a dipole

# Arguments
- `k`: wave number
- `coords`: dipole coordinates (N x 3 matrix)
- `P`: polarizations (vector of length 3N; Px1,Py1,Pz1 ... PxN,PyN,PzN)
- `r`: coord for the point at which to calculate the near-field
"""
function E_sca(k, coords, P, r)

    E_sca = zero(eltype(P))
    for i in eachindex(coords)
        A = calc_Ajk(k, coords[i], r)
        E_sca += A * P[:, i]
    end
    return E_sca
end


@doc raw"""
    E_sca_far_field(k, coords, P, r)

In the far-field, the scattered electric field is given by:
```math
E_{sca} = \frac{k^2 \exp(ikr)}{r} \sum \limits_{j=1}^{N} \exp(-ik\hat{r} \cdot r_j)(\hat{r}\hat{r}-1_3)P_j
```
Note: coordinates are relative to origin

# Arguments
- `k`: wave number
- `coords`: dipole coordinates (N x 3 matrix)
- `P`: polarizations (vector of length 3N; Px1,Py1,Pz1 ... PxN,PyN,PzN)
- `r`: coord for the point at which to calculate the far-field

"""
function E_sca_far_field(k, coords, P, r)

    r_norm = norm(r)
    r_hat = r / r_norm

    k^2 * exp(im * k * r_norm) / r_norm * sum(i ∈ eachindex(coords)) do
        exp(-im * k * dot(r_hat, coords[i])) * (r_hat' .* r_hat - I) * P[:, i]
    end
end


# function E_sca_far_field(k, coords, P, r)
#
#     r_norm = norm(r)
#     r_hat = r / r_norm
#
#     E = zeros(eltype(P))
#     for i in eachindex(coords)
#         E += exp(-im * k * dot(r_hat, coords[i])) * (r_hat' .* r_hat - I) * P[:,i]
#     end
#
#     return  k^2 * exp(im * k * r_norm) / r_norm * E
#
#     # return k^2 * exp(im * k * r_norm) / r_norm *
#     # sum(exp(-im * k * dot(r_hat, coords[i])) * (r_hat' .* r_hat - I) * P[:, i] for i ∈ eachindex(coords))
#
# end
