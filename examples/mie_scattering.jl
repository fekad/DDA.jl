# using Plots:append!
using SpecialFunctions
using Plots
using LaTeXStrings

plotlyjs()


@doc raw"""
    rayleigh_scattering(m, wavelength, diameter, n_medium=1.0)

Computes Mie efficencies of a spherical particle in the Rayleigh regime (``x = \pi \, d_p / \lambda \ll 1``) given refractive index `m`, `wavelength`, and `diameter`. Uses Rayleigh-regime approximations:
```math
Q_{sca}=\frac{8x^4}{3}\left|{\frac{m^2-1}{m^2+2}}\right|^2}
```
```math
Q_{abs}=4x\:\text{Im}\left\{\frac{m^2-1}{m^2+2}\right\}}
```
```math
Q_{ext}=Q_{sca}+Q_{abs}}
```
```math
Q_{back}=\frac{3Q_{sca}}{2}}
```
```math
Q_{ratio}=1.5}
```
```math
Q_{pr}=Q_{ext}}
```

# Arguments:
- m: The complex refractive index, with the convention ``m = n+ik``.
- wavelength: The wavelength of incident light, in nanometers.
- diameter: The diameter of the particle, in nanometers.
- n_medium: The refractive index of the surrounding medium. This must be positive, nonzero, and real. Any imaginary part will be discarded.

# Returns
- qext, qsca, qabs, qback: The Mie efficencies described above.

# Examples
```julia-repl
julia> rayleigh_scattering(1.33+0.01im, 870., 50.)

(0.0041753430994240295, 0.00011805645915412197, 0.004057286640269908, 0.00017708468873118297)
```
"""
function rayleigh_scattering(m::ComplexF64, wavelength::Float64, diameter::Float64, n_medium::Float64=1.0)
    m /= n_medium
    wavelength /= n_medium

    x = pi * diameter / wavelength

    if x ≈ 0
        return  0, 0, 0, 0
    end

    LL = (m^2 - 1) / (m^2 + 2) # Lorentz-Lorenz term

    qsca = 8 * x^4 / 3 * abs(LL)^2  # B&H eq 5.8
    qabs = 4 * x * imag(LL) # B&H eq. 5.11
    qext = qsca + qabs
    qback = 3 * qsca / 2 # B&H eq. 5.9

    # qratio = 1.5
    # qpr = qext

    return qext, qsca, qabs, qback
end


rayleigh_scattering(1.33 + 0.01im, 870., 0.)
rayleigh_scattering(1.33 + 0.01im, 870., 50.)

lambda = LinRange(600, 1000, 1001)
qext =  [rayleigh_scattering(1.33 + 0.01im, l, 50.)[1] for l in lambda]
plot(lambda, qext)



@doc raw"""
    mie_scattering(m, wavelength, diameter, n_medium=1.0)

Computes Mie efficencies *Q* and asymmetry parameter *g* of a single, homogeneous particle. Uses `mie_ab` to calculate ``a_n`` and ``b_n``, and then calculates *Q* via:
```math
Q_{ext}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)\:\text{Re}\left\{a_n+b_n\right\}}
```
```math
Q_{sca}=\frac{2}{x^2}\sum_{n=1}^{n_{max}}(2n+1)(|a_n|^2+|b_n|^2)}
```
```math
Q_{abs}=Q_{ext}-Q_{sca}}
```
```math
Q_{back}=\frac{1}{x^2} \left| \sum_{n=1}^{n_{max}}(2n+1)(-1)^n(a_n-b_n) \right| ^2}
```
```math
Q_{ratio}=\frac{Q_{back}}{Q_{sca}}}
```
```math
g=\frac{4}{Q_{sca}x^2}\left[\sum\limits_{n=1}^{n_{max}}\frac{n(n+2)}{n+1}\text{Re}\left\{a_n a_{n+1}^*+b_n b_{n+1}^*\right\}+\sum\limits_{n=1}^{n_{max}}\frac{2n+1}{n(n+1)}\text{Re}\left\{a_n b_n^*\right\}\right]}
```
```math
Q_{pr}=Q_{ext}-gQ_{sca}}
```
where asterisks denote the complex conjugates.

# Arguments
- m: The complex refractive index, with the convention *m = n+ik*.
- wavelength: The wavelength of incident light, in nanometers.
- diameter: The diameter of the particle, in nanometers.
- n_medium: The refractive index of the surrounding medium. This must be positive, nonzero, and real. Any imaginary part will be discarded.

# Returns
- qext, qsca, qabs, g, qpr, qback, qratio: The Mie efficencies described above.

# Examples

For example, compute the Mie efficencies of a particle ``300 nm`` in diameter with ``m = 1.77+0.63i``, illuminated by ``\lambda = 375 nm``:
```julia-repl
julia> mie_scattering(1.77+0.63im, 375, 300)

(2.858497199156411, 1.3149276685170943, 1.543569530639317, 0.2014551048135256)
```
"""
function mie_scattering(m, wavelength, diameter, n_medium=1.0)
    m /= n_medium
    wavelength /= n_medium

    x = pi * diameter / wavelength

    if x ≈ 0
        return 0., 0., 0., 0.
    end

    nmax = round(Int, 2 + x + 4 * (x^(1 / 3)))

    an, bn = mie_ab(m, x, nmax)

    qext = 2 / x^2 * sum((2n + 1) * real(an[n] + bn[n]) for n = 1:nmax)
    qsca = 2 / x^2 * sum((2n + 1) * (abs(an[n])^2 + abs(bn[n])^2) for n = 1:nmax)
    qabs = qext - qsca
    qback = 1 / x^2 * abs(sum((2n + 1) * (-1)^n * (an[n] - bn[n]) for n = 1:nmax))^2

    # g = 4 / (qsca * x^2) * (sum(n * (n + 2) / (n + 1) * real(an[n]*an[n+1]'+bn[n]*bn[n+1]') for n = 1:nmax-1) + sum((2n + 1) / (n * (n + 1)) * real(an[n]*bn[n]')  for n = 1:nmax-1))
    # qpr = qext - qsca * g
    # qratio = qback / qsca

    return qext, qsca, qabs, qback
end

@doc raw"""
    mie_ab(m, x, nmax)

Computes external field coefficients ``a_n`` and ``b_n`` based on inputs of `m` and ``x = \pi \, d_p / \lambda``.


# Arguments
- m: The complex refractive index with the convention ``m = n+ik``.
- x: The size parameter ``x = \pi \, d_p / \lambda``.

# Returns
- an, bn: Arrays of size ``n_max = 2 + x + 4x^{1/3}``

"""
function mie_ab(m, x, nmax)
    mx = m * x
    nmx = max(nmax, round(Int, abs(mx))) + 16

    n = 1:nmax
    nu = n .+ 0.5
    sx = sqrt(0.5 * pi * x)

    px = sx * besselj.(nu, x)


    p1x = append!([sin(x)], px[1:end - 1])

    chx = -sx * bessely.(nu, x)
    ch1x = append!([cos(x)], chx[1:end - 1])

    gsx = px - (0 + 1im) * chx #
    gs1x = p1x - (0 + 1im) * ch1x #

    # B&H Equation 4.89
    Dn = zeros(ComplexF64, nmx)
    for i in nmx - 1:-1:1
        Dn[i] = (i / mx) - (1 / (Dn[i + 1] + i / mx))
    end

    D = Dn[2:nmax + 1] # Dn(mx), drop terms beyond nMax

    da = D / m + n / x
    db = m * D + n / x

    an = @. (da * px - p1x) / (da * gsx - gs1x)
    bn = @. (db * px - p1x) / (db * gsx - gs1x)

    return an, bn
end


mie_scattering(1.33 + 0.1im, 870., 0.)
mie_scattering(1.33 + 0.1im, 870., 50.)
mie_scattering(1.5 + 0.5im,532,200)


lambda = LinRange(10, 1000, 1001)
qext =  [mie_scattering(1.33 + 0.0im, l, 50.)[1] for l in lambda]
plot(lambda, qext)



# m1 = 1.33; # relative refractive index of water
m1 = 1.33 + .1i;   # imag. component to demonstrate absorption
k = 2 * pi;          # wave number
d = 1 / (abs(m1) * k); # lattice spacing





@doc raw"""
    scattering_function(m, wavelength, diameter, n_medium=1.0, minAngle=0, maxAngle=180, angularResolution=0.5)

Creates arrays for plotting the angular scattering intensity functions in theta-space with parallel, perpendicular, and unpolarized light. Also includes an array of the angles for each step. This angle can be in either degrees, radians, or gradians for some reason. The angles can either be geometrical angle or the qR vector (see `Sorensen, M. Q-space analysis of scattering by particles: a review. J. Quant. Spectrosc. Radiat. Transfer 2013, 131, 3-12 <http://www.sciencedirect.com/science/article/pii/S0022407313000083>`_). Uses `mie_S1_S2` to compute ``S_1`` and ``S_2``, then computes parallel, perpendicular, and unpolarized intensities by

```math
SL(\theta)=|S_1|^2}
```
```math
SR(\theta)=|S_2|^2}
```
```math
SU(\theta)=\frac{1}{2}(SR + SL)}
```

# Arguments

- theta: An array of the angles used in calculations.
- m: The complex refractive index with the convention *m = n+ik*.
- wavelength: The wavelength of incident light, in nanometers.
- diameter: The diameter of the particle, in nanometers.
- n_medium: The refractive index of the surrounding medium. This must be positive, nonzero, and real. Any imaginary part will be discarded.

# Returns

- SL: An array of the scattered intensity of left-polarized (perpendicular) light. Same size as the **theta** array.
- SR: An array of the scattered intensity of right-polarized (parallel) light. Same size as the **theta** array.
- SU: An array of the scattered intensity of unpolarized light, which is the average of SL and SR. Same size as the **theta** array.
"""
function scattering_function(theta, m, wavelength, diameter, n_medium=1.0)

    m /= n_medium
    wavelength /= n_medium

    x = pi * diameter / wavelength

    n = length(theta)

    if x ≈ 0
        return 0., 0., 0.
    end

    SL = zeros(n)
    SR = zeros(n)
    SU = zeros(n)
    for i in 1:n
        mu = cos(theta[i])
        S1, S2 = mie_S1_S2(m, x, mu)
        SL[i] = real(sum(S1' * S1))
        SR[i] = real(sum(S2' * S2))
        SU[i] = (SR[i] + SL[i]) / 2
    end

    return SL, SR, SU
end

@doc raw"""
    mie_S1_S2(m,x,mu)

Calculates ``S_1`` and ``S_2`` at ``\mu = cos(\theta)``, where ``\theta`` is the scattering angle.

Uses `mie_ab` to calculate `a_n` and `b_n`, and `mie_pi_tau` to calculate ``\pi_n`` and ``\tau_n``. ``S_1`` and ``S_2`` are calculated by:

```math
S_1 = \sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)} (a_n \pi_n + b_n \tau_n)
```
```math
S_2 = \sum\limits_{n=1}^{n_{max}} \frac{2n+1}{n(n+1)} (a_n \tau_n + b_n \pi_n)
```

# Parameters
- m: The complex refractive index with the convention *m = n+ik*.
- x: The size parameter ``x = \pi \, d_p / \lambda``.
- mu: The cosine of the scattering angle.

# Returns
-S1, S2: The ``S_1`` and ``S_2`` values.

"""
function mie_S1_S2(m, x, mu)

    nmax = round(Int, 2 + x + 4 * x^(1 / 3))

    an, bn = mie_ab(m, x, nmax)
    pin, taun = mie_pi_tau(mu, nmax)

    S1 = sum((2n + 1) / (n * (n + 1)) * (an[n] * pin[n] + bn[n] * taun[n]) for n = 1:nmax)
    S2 = sum((2n + 1) / (n * (n + 1)) * (an[n] * taun[n] + bn[n] * pin[n]) for n = 1:nmax)

    return S1, S2
end


@doc raw"""
    mie_pi_tau(mu, nmax)

Calculates ``\pi_n`` and ``\tau_n``.

This function uses recurrence relations to calculate ``\pi_n`` and ``\tau_n``, beginning with ``\pi_0 = 1``, ``\pi_1 = 3 \mu`` (where ``mu`` is the cosine of the scattering angle), ``\tau_0 = \mu``, and ``\tau_1 = 3 cos(2 cos^{-1}(\mu))``:
```math
\pi_n = \frac{2n-1}{n-1} \mu \pi_{n-1} - \frac{n}{n-1} \pi_{n-2}
```
```math
\tau_n = n \mu\ pi_n - (n+1) \ pi_{n-1}
```

# Arguments
- mu : The cosine of the scattering angle.
- nmax: The number of elements to compute. Typically, ``n_{max} = floor(2+x+4x^{1/3})``, but can be given any integer.

# Returns
- p, t: The ``\pi_n`` and ``\tau_n`` arrays, of length `nmax`.
"""
function mie_pi_tau(mu, nmax)

    p = zeros(nmax)
    t = zeros(nmax)

    p[1] = 1
    p[2] = 3mu

    t[1] = mu
    t[2] = 3 * cos(2 * acos(mu))

    for n in 3:nmax
        p[n] = (2n - 1) / (n - 1) * mu * p[n - 1] - n / (n - 1) * p[n - 2]
        t[n] = n * mu * p[n] - (n + 1) * p[n - 1]
    end

    return p, t
end




m = 1.7 + 0.5im
w = 532
d = 5000
theta = LinRange(0, pi, 361)
SL, SR, SU = scattering_function(theta, m, w, d)

plot(
    theta, SL,
    label="Parallel Polarization",
    yaxis=:log,
    ylabel="Intensity (\$|S|^2}\$)",
    xlabel="ϴ",
    xticks=([0, pi / 4, pi / 2, 3 * pi / 4, pi],
        ["0", L"\frac{\pi}{4}", L"\frac{\pi}{2}", L"\frac{3\pi}{4}", L"\pi"]
    )
)
plot!(theta, SR, yaxis=:log, label="Perpendicular Polarization")
plot!(theta, SU, yaxis=:log, label="Unpolarized")





m1 = 1.33; # relative refractive index of water
# m1 = 1.33 + .1im   # imag. component to demonstrate absorption
k = 2 * pi          # wave number
d = 1 / (abs(m1) * k) # lattice spacing

# number of dipoles in the approximate sphere; more is required as the
# radius increases
nrange = [8, 32, 136, 280, 552, 912, 1472, 2176, 3112, 4224, 5616, 7208, 9328, 11536];
# nrange = [8, 32, 136, 280, 552, 912, 1472]

# the corresponding effective radii of the spheres
a_eff = (3 * nrange / 4π).^(1 / 3) * d

arange = range(a_eff[1], a_eff[end], length=1001)

λ = 2π / k
Q = zeros(length(arange), 3)
for (i, a) = enumerate(arange)
    # @show (i, a)
    r = mie_scattering(m1, λ, 2.1 * a)
    Q[i,:] .= r[1:3]
end

# Here, we plot the efficiencies Q instead of the cross sections C
# Q = C/(pi*r^2)
plot(k * arange, Q;
    labels=["Q_{ext}" "Q_{abs}" "Q_{scat}"],
    title="m = $m1"
    )
ylabel!("Q")
xlabel!("2 π a λ") # size parameter ka



nothing


# m1 = 1.33; # relative refractive index of water
m1 = 1.33 + .1im   # imag. component to demonstrate absorption
k = 2 * pi          # wave number
d = 1 / (abs(m1) * k) # lattice spacing

# number of dipoles in the approximate sphere; more is required as the
# radius increases
nrange = [8, 32, 136, 280, 552, 912, 1472, 2176, 3112, 4224, 5616, 7208, 9328, 11536];
# nrange = [8, 32, 136, 280, 552, 912, 1472]

# the corresponding effective radii of the spheres
a_eff = (3 * nrange / 4π).^(1 / 3) * d

arange = range(a_eff[1], a_eff[end], length=1001)

λ = 2π / k
Q = zeros(length(arange), 3)
for (i, a) = enumerate(arange)
    # @show (i, a)
    r = mie_scattering(m1, λ, 2.1 * a)
    Q[i,:] .= r[1:3]
end

# Here, we plot the efficiencies Q instead of the cross sections C
# Q = C/(pi*r^2)
plot(k * arange, Q;
    labels=["Q_{ext}" "Q_{abs}" "Q_{scat}"],
    title="m = $m1"
    )
ylabel!("Q")
xlabel!("2 π a λ") # size parameter ka








