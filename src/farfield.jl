# TODO:single alpha
@doc raw"""
The absorption cross section:
$$
C_{abs} = \frac{4 \pi k}{|E_0|^2} \sum \limits_{j=1}^{N} \left\{ \mathrm{Im}(P_j \cdot (\alph_j^{-1}) \cdot P_j^*) - \frac{2}{3} k^3|P_j|^2 \right\}
$$
"""
function C_abs(k, E0, P, alph)
    c = zero(Float64)
    for i in 1:length(P)
        c += -imag(dot(P[i], 1 / alph, P[i])) - 2 / 3 * k^3 * norm(P[i])^2
    end
    return 4π * k / norm(E0)^2 * c
end

C_abs(sol::GridSolution) = C_abs(norm(sol.prob.Einc.kvec), norm(sol.prob.Einc.E₀), sol.P, sol.alphas)



@doc raw"""
The extinction cross section:
$$
C_{ext} = \frac{4 \pi k}{|E_0|^2} \sum \limits_{j=1}^{N} \mathrm{Im}(\vec{E}_{inc, j}^* \cdot P_j)
$$
"""
# C_ext(k, E0, Ei, P) = 4π * k / norm(E0)^2 * imag(dot(Ei, P))
function C_ext(k, E0, Einc, P)
    c = zero(Float64)
    for j = 1:length(P)
        c += imag(dot(Einc[j], P[j]))
    end
    return 4π * k / norm(E0)^2 * c
end

C_ext(sol::GridSolution)= C_ext(norm(sol.prob.Einc.kvec), norm(sol.prob.Einc.E₀), sol.Einc, sol.P)

@doc raw"""
The scattering cross section
$$
C_{sca} = C_{ext} - C_{abs}
$$
"""
C_sca(k, E0, Einc, P, alph) = C_ext(k, E0, Einc, P) - C_abs(k, E0, P, alph)

C_sca(sol::GridSolution) = C_ext(sol) - C_abs(sol)
