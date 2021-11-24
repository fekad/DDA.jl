@doc raw"""
The absorption cross section:
$$
C_{abs} = \frac{4 \pi k}{|E_0|^2} \sum \limits_{j=1}^{N} \left\{ \mathrm{Im}(P_j \cdot (\alph_j^{-1}) \cdot P_j^*) - \frac{2}{3} k^3|P_j|^2 \right\}
$$
"""
function C_abs(k, E0, Ei, P, alph)
    A = 4π * k / norm(x)^2
    N = size(alph,2)
    # @show A

    C = 0
    for i in 1:N
        C += -imag(dot(P[:,i], 1 / alph[i], P[:,i])) - 2 / 3 * k^3 * norm(P[:,i])^2
        # C += 4π * k / norm(E0)^2 * sum(-imag(P[j]*(1/alph[j]*I)*P[j]') - 2/3 *k^3 * P[j]^2 )
        # !!!!!!!!!!!!!!
        # !dot conjugate ???
    end
    return A * C
end

@doc raw"""
The extinction cross section:
$$
C_{ext} = \frac{4 \pi k}{|E_0|^2} \sum \limits_{j=1}^{N} \mathrm{Im}(\vec{E}_{inc, j}^* \cdot P_j)
$$
"""
C_ext(k, E0, Ei, P) = 4π * k / norm(E0)^2 * imag(dot(Ei, P))

@doc raw"""
The scattering cross section
$$
C_{sca} = C_{ext} - C_{abs}
$$
"""
C_sca(k, E0, Ei, P, alph) = C_ext(k, E0, Ei, P) - C_abs(k, E0, Ei, P, alph)
