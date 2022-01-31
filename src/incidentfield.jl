
# function field(::AbstractIncidentField, r::V) where {Dim,V<:AbstractVector{<:StaticVector{Dim,<:Real}}}
# function field(::AbstractIncidentField, r::AbstractVector{<:StaticVector{Dim,T}) where {Dim,T<:Real}}

function field(f::AbstractIncidentField, r::AbstractVector{<:AbstractVector{<:Real}})
    E = similar(r, SVector{3,Complex{Float64}})
    for i = 1:length(r)
        E[i] = field(f, r[i])
    end
    E
end


raw"""
# Incident electric field

Incident planewave:
```math
E_{inc}(r) = E_0 e^{i (kr - \omega t)}
```

Arbitrary incident electric field:
- [Jones vector](https://en.m.wikipedia.org/wiki/Jones_calculus) in the lab frame (x-y polarisation plane, z propagation),
- 3D rotation of the coordinate system

Visualisation:
- Polarisation state can be represented as a point on the Bloch (Poincare) sphere.
"""


@doc raw"""
    PlaneWave(k, e, θ, ϕ)

Suppose that a monochromatic plane wave of light is travelling in the positive ``z``-direction, with angular frequency ``\omega`` and wave vector ``\mathbf{k} = (0,0,k)``, where the wavenumber k = ω/c. Then the electric and magnetic fields E and H are orthogonal to k at each point; they both lie in the plane "transverse" to the direction of motion. Furthermore, H is determined from E by 90-degree rotation and a fixed multiplier depending on the wave impedance of the medium. So the polarization of the light can be determined by studying E. The complex amplitude of E is written

```math
\begin{pmatrix} E_{x}(t) \\ E_{y}(t) \\ 0 \end{pmatrix} =
\begin{pmatrix} E_{0x} e^{i(kz - \omega t + \phi_{x})} \\ E_{0y} e^{i(kz - \omega t + \phi_{y})} \\ 0 \end{pmatrix} =
\begin{pmatrix} E_{0x} e^{i\phi_{x}} \\ E_{0y} e^{i\phi_{y}} \\ 0 \end{pmatrix} e^{i(kz - \omega t)}.
```

Linear polarized in the x direction (typically called "horizontal"):
```math
|H\rangle = \begin{pmatrix}1\\0\end{pmatrix}
```

Linear polarized in the y direction (typically called "vertical"):
```math
|V\rangle = \begin{pmatrix}0\\1\end{pmatrix}
```

Linear polarized at 45° from the x axis (typically called "diagonal" L+45):
```math
|D\rangle = \frac{1}{\sqrt{2}} \left(|H\rangle +|V\rangle \right) = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\1\end{pmatrix}
```

Linear polarized at -45° from the x axis (typically called "anti-diagonal" L-45):
```math
|A\rangle =\frac{1}{\sqrt{2}} \left(|H\rangle - |V\rangle\right) = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\-1\end{pmatrix}
```

Right-hand circular polarized (typically called "RCP" or "RHCP"):
```math
|R\rangle =\frac{1}{\sqrt{2}} \left(|H\rangle - i|V\rangle\right) = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\-i\end{pmatrix}
```

Left-hand circular polarized (typically called "LCP" or "LHCP"):
```math
|L\rangle =\frac{1}{\sqrt{2}} \left(|H\rangle + i|V\rangle\right) = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\+i\end{pmatrix}
```

Arbitrary direction:
```math
\mathbf{E}(\mathbf{r}, t) =
\mathbf{R} \begin{pmatrix} E_{0x} e^{i\phi_{x}} \\ E_{0y} e^{i\phi_{y}} \\ 0 \end{pmatrix} e^{i (\mathbf{R} \mathbf{k}) \cdot \mathbf{r} }e^{-i \omega t}.
```
"""
struct PlaneWave <: AbstractIncidentField
    kvec::SVector{3,<:Real}
    E₀::SVector{3,Complex{<:Real}}
end

# "minimal" representation of a plane wave
function PlaneWave(k::Real, e::StaticVector{2,<:Complex{<:Real}}, θ::Real, ϕ::Real)
    R = RotZY(θ, ϕ)
    E₀ = R[:, 1:2] * e
    kvec = R[:, 3] * k
    return PlaneWave(kvec, E₀)
end

PlaneWave(k, e::AbstractArray, θ, ϕ) = PlaneWave(k, convert(SVector{2,ComplexF64}, e), θ, ϕ)

wavenumber(f::PlaneWave) = norm(f.kvec)
wavevector(f::PlaneWave) = f.kvec
E0(f::PlaneWave) = norm(f.E₀)
polarisation(f::PlaneWave) = normalize(f.E₀)

field(f::PlaneWave, r::AbstractVector{<:Real}) = f.E₀ * exp(im * dot(f.kvec, r))



struct EELS <: AbstractIncidentField
    k::Float64 # wavenumber
    v::Float64 # charge velocity (propagate along z direction)
    origin::SVector{2,Float64} # position in the xy plane
end

wavenumber(f::EELS) = f.k
wavevector(f::EELS) = SVector{3,Float64}(0, 0, f.k)

K₀(z) = besselk(0, z)
K₁(z) = besselk(1, z)

function field(f::EELS, r::AbstractVector{<:Real})
    d = r[1:2] - f.origin
    dnorm = norm(d)

    A = q.val / (2π * ε_0.val) * f.k / f.v * exp(1im * f.k * r[3])

    return [
        A * d[1] / dnorm * K₁(f.k * dnorm),
        A * d[2] / dnorm * K₁(f.k * dnorm),
        -A * 1im * K₀(f.k * dnorm)
    ]
end

