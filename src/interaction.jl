
@doc raw"""
    calc_Ajk(k, r_jk)

The electric field from a radiating electric dipole:
```math
E = \frac{1}{4 \pi \varepsilon_0} \left\{
    k^2 (\hat{r} \times p) \times \hat{r} \frac{e^{ikr}}{r} +
    \left[3\hat{r}(\hat{r}\cdot p)-p\right] \left( \frac{1}{r^3}-\frac{ik}{r^2} \right) e^{ikr}
\right\}
```

`A_{jk}` is the tensor that represents the interaction between a receiving dipole at `r_j` and the radiating dipole at `r_k`.
```math
A_{jk} = \frac{\exp(i k r_{jk})}{r_{jk}} \left[
    k^2 (\hat{r}_{jk} \hat{r}_{jk} - 11) +
    \frac{i k r_{jk} - 1}{r_{jk}^2} (3 \hat{r}_{jk} \hat{r}_{jk} - 11)
\right]
```
"""
function calc_Ajk(k, r_jk)
    r = norm(r_jk)
    rn = r_jk / r

    A_jk = exp(im * k * r) / r * (k^2 * (rn * rn' - I) + (im * k * r - 1) / r^2 * (3 * rn * rn' - I))
    return A_jk
end

calc_Ajk(k,  r_j, r_k) = calc_Ajk(k, r_j - r_k)

