@doc raw"""
    C_abs(k, E0, P, alpha)

The absorption cross section:
```math
C_{abs} = \frac{4 \pi k}{|E_0|^2} \sum \limits_{j=1}^{N} \left\{ \mathrm{Im}(P_j \cdot (\alpha_j^{-1}) \cdot P_j^*) - \frac{2}{3} k^3|P_j|^2 \right\}
```
"""
function C_abs(k, E0, P, alpha)
    T = Float64
    c = zero(T)
    for i in eachindex(P)
        c += -imag(dot(P[i], inv(alpha[i]), P[i])) - 2 / 3 * k^3 * norm(P[i])^2
    end
    return 4π * norm(k) / norm(E0)^2 * c
end


@doc raw"""
    C_ext(k, E0, Einc, P)

The extinction cross section:
```math
C_{ext} = \frac{4 \pi k}{|E_0|^2} \sum \limits_{j=1}^{N} \mathrm{Im}(\vec{E}_{inc, j}^* \cdot P_j)
```
"""
function C_ext(k, E0, Einc, P)
    T = Float64
    c = zero(T)
    for j in eachindex(P)
        c += imag(dot(Einc[j], P[j]))
    end
    return 4π * norm(k) / norm(E0)^2 * c
end


@doc raw"""
    C_sca(k, E0, Einc, P, alpha)

The scattering cross section:
```math
C_{sca} = C_{ext} - C_{abs}
```
"""
C_sca(k, E0, Einc, P, alpha) = C_ext(k, E0, Einc, P) - C_abs(k, E0, P, alpha)