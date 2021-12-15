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


abstract type IncindentField <: AbstractField end


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

Linear polarized at −45° from the x axis (typically called "anti-diagonal" L−45):
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
struct PlaneWave <: IncindentField
    k::Float64
    e::SVector{2,ComplexF64}
    θ::Float64
    ϕ::Float64
end

# E = PlaneWave(1, (1, 1im), 1, 1)

function params(E::PlaneWave)
    R = RotZY(E.θ, E.ϕ)
    E₀ = R[:,1:2] * E.e
    kvec = R[:,3] * E.k
    return kvec, E₀
end


function field(E::PlaneWave, r::AbstractVector)
    return E₀ * exp(im * dot(kvec, r))
end


struct PlaneWave2 <: IncindentField
    kvec::SVector{3,Float64}
    E₀::SVector{3,ComplexF64}
end

# syntax sugar
function PlaneWave2(k ,e, θ, ϕ)
    R = RotZY(θ, ϕ)
    E₀ = R[:,1:2] * e
    kvec = R[:,3] * k
    return PlaneWave2(kvec, E₀)
end

function field(f::PlaneWave2, r::AbstractVector{Point3{T}}) where T
    E = similar(r, SVector{3, ComplexF64})
    for i in 1:length(r)
        E[i] = field(f, r[i])
    end
    E
end


field(f::PlaneWave2, r::Point3) = f.E₀ * exp(im * dot(f.kvec, r))

# function field(E::PlaneWave, r::AbstractVector{AbstractVector})
#     out = Array{ComplexF64}(undef, (3, length(r)))
#     for i in eachindex(r)
#         out[:,i] = field(E, r[i])
#     end
#     return out
# end


function field(E::PlaneWave, r::AbstractArray)

    # @boundcheck if size(r, 1) != 3
    #     throw("dimension error")
    # end

    out = Array{ComplexF64}(undef, (3, length(r)))
    for i in eachindex(r)
        out[:,i] = field(E, r[i])
    end
    return out
end





# """
# Incident wave:
# ```math
# E_{inc,j} = E_0 \exp(-i k r_j - i \omega t)
# ```
# """
#
# https://discourse.julialang.org/t/struct-of-arrays-soa-vs-array-of-structs-aos/30015/16
#
# # polarization:
# # left and right hand circular polarized VS the orthonormal Cartesian vectors ex and ey
# struct PlaneWave <: Field
#     E₀::SVector{3,Float64}
#     k::SVector{3,Float64}
# end
#
# function field(r, k, E₀)
#     return E₀ * exp(-im * dot(k, r))
# end
#
# function field(r, E::PlaneWave)
#     return field(r, E.E₀, E.k)
# end
#
# function field(r::Vector{SVector{3, Float64}}, e::PlaneWave)
#     out = Array{ComplexF64}(undef, (3, length(r)))
#     for i in eachindex(r)
#         out[:,i] = field(r[i], e)
#     end
#     return out
# end
#
# # E_inc(r, E::PlaneWave) = E.E0 * exp(-im * E.k)
# # E_inc_td(t, r, E::PlaneWave, omega) = field(E) * exp(im * omega * t)
# #
# # k = [0 0 1]
# # # u = [0,0,1] polarisibility
# # E_inc(r, E_0, k) = E_0 * exp(-im * dot(k,r))
# #
#
#
# # function E_inc(E0, kvec, r)
# #     Ei = zeros(ComplexF64, 3, length(r))
# #
# #     for (i, ri) in enumerate(r)
# #         Ei[:,i] = E0 .* exp.(im * dot(kvec, ri))
# #     end
# #
# #     return reshape(Ei, :)
# # end




