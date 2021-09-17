
@doc raw"""
The electric field from a radiating electric dipole:
$$
E = \frac{1}{4 \pi \varepsilon_0} \left\{
    k^2(\hat{r} \times p) \times \hat{r} \frac{e^{ikr}}{r} +
    \left[3\hat{r}(\hat{r}\cdot p)-p\right] \left( \frac{1}{r^3}-\frac{ik}{r^2} \right) e^{ikr}
\right\}
$$

$A_{jk}$ is the tensor that represents the interaction between a receiving dipole at $r_j$ and the radiating dipole
at $r_k$.
$$
A_{jk} = \frac{\exp(i k r_{jk})}{r_{jk}} \left[
    k^2(\hat{r}_{jk} \hat{r}_{jk} - 11) + \frac{i k r_{jk} - 1}{r_{jk}^2} (3 \hat{r}_{jk} \hat{r}_{jk} - 11)
\right]
$$
"""
function calc_Ajk(k, rj, rk)
    r = rj - rk
    r_jk = norm(r)
    r_n = r / r_jk

    A_jk = exp(im * k * r_jk) / r_jk * (k^2 * (r_n * r_n' - I) + (im * k * r_jk - 1) / (r_jk^2) * (3 * r_n * r_n' - I))
    return A_jk
end

calc_Ajj(alpha) = Diagonal(1 / alpha * I, 3)



function interaction(k, r, alpha)
    N = length(r)
    out = zeros(ComplexF64, 3, N, 3, N)

    for i in 1:N
        # out[i,i,:,:] = 1/alpha[i] * I
        out[:,i,:,i] = DDA.calc_Ajj(alpha[i])
    end

    for i in 2:N
        for j in 1:i - 1
            # @show i, j
            # @show k, r[i], r[j]
            # @show DDA.calc_Ajk(k, r[i], r[j])
            out[:,i,:,j] = DDA.calc_Ajk(k, r[i], r[j])
        end
    end

    # DDA.calc_Ajk(k, r[2], r[1])
    return Symmetric(reshape(out, 3N, 3N), :L)
end

@doc raw"""
The absorption cross section:
$$
C_{abs} = \frac{4 \pi k}{|E_0|^2} \sum \limits_{j=1}^{N} \left\{ \mathrm{Im}(P_j \cdot (\alpha_j^{-1}) \cdot P_j^*) - \frac{2}{3} k^3|P_j|^2 \right\}
$$
"""
function C_abs(k, E0, Ei, P, alpha)
    A = 4π * k / norm(E0)^2
    N = length(alpha)
    @show A

    C = 0
    for i in 1:N
        C += -imag(dot(P[:,i], 1 / alpha[i], P[:,i])) - 2 / 3 * k^3 * norm(P[:,i])^2
        # C += 4π * k / norm(E0)^2 * sum(-imag(P[j]*(1/alpha[j]*I)*P[j]') - 2/3 *k^3 * P[j]^2 )
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
C_sca(k, E0, Ei, P, alpha) = C_ext(k, E0, Ei, P) - C_abs(k, E0, Ei, P, alpha)


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





