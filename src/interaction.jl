
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
function calc_Ajk(k, r_jk)
    r = norm(r_jk)
    rn = r_jk / r

    A_jk = exp(im * k * r) / r * (k^2 * (rn * rn' - I) + (im * k * r - 1) / r^2 * (3 * rn * rn' - I))
    return A_jk
end

function calc_Ajk(k, rj, rk)
    r_jk = rj - rk
    return calc_Ajk(k, r_jk)
end

# TODO: allpha is a 3x3 (symmetric) tensor by defualt )
calc_Ajj(alph) = Diagonal(1 / alph * I, 3)



function interactions(k, r, alph)
    N = length(r)
    out = zeros(ComplexF64, 3, N, 3, N)

    for i in 1:N
        # out[i,i,:,:] = 1/alph[i] * I
        out[:,i,:,i] = calc_Ajj(alph[i])
    end

    for i in 2:N
        for j in 1:i - 1
            # @show i, j
            # @show k, r[i], r[j]
            # @show DDA.calc_Ajk(k, r[i], r[j])
            out[:,i,:,j] = calc_Ajk(k, r[i], r[j])
        end
    end

    # DDA.calc_Ajk(k, r[2], r[1])
    return Symmetric(reshape(out, 3N, 3N), :L)
end






